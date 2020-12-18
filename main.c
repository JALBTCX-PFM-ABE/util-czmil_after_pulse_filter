
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

 /********************************************************************
 *
 * Module Name : main.c
 *
 * Author/Date : Jan C. Depner - 10/27/16
 *
 * Description : Filters out after pulse near strongest return
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <math.h>


/* Local Includes. */

#include "nvutility.h"

#include "czmil.h"

#include "version.h"


void usage ()
{
  fprintf (stderr, "\nUsage: czmil_after_pulse_filter [-L LNS] [-U UNS] [-d DIFF] [-1] [-2] [-3] [-4] [-5] [-6] [-7] [-9] CPFFILE\n");
  fprintf (stderr, "Where:\n");
  fprintf (stderr, "\tLNS = Nanoseconds from strongest return to after pulse, lower bound (20)\n");
  fprintf (stderr, "\tUNS = Nanoseconds from strongest return to after pulse, upper bound (40)\n");
  fprintf (stderr, "\tDIFF = Minimum amplitude difference between the strongest return interest point\n");
  fprintf (stderr, "\t\tand the after-pulse return interest point.  If this threshold is exceeded and\n");
  fprintf (stderr, "\t\tthe after-pulse interest point is between the LNS and UNS the point will be\n");
  fprintf (stderr, "\t\tinvalidated.\n\n");
  fprintf (stderr, "\t-1 = filter channel 1\n");
  fprintf (stderr, "\t-2 = filter channel 2\n");
  fprintf (stderr, "\t-3 = filter channel 3\n");
  fprintf (stderr, "\t-4 = filter channel 4\n");
  fprintf (stderr, "\t-5 = filter channel 5\n");
  fprintf (stderr, "\t-6 = filter channel 6\n");
  fprintf (stderr, "\t-7 = filter channel 7\n");
  fprintf (stderr, "\t-8 = filter channel 9\n\n");
  fprintf (stderr, "\tThis filter only works for LAND processed data.\n\n");
  exit (-1);
}



int32_t main (int32_t argc, char **argv)
{
  char               cpf_file[1024], cwf_file[1024];
  int32_t            cpf_hnd = -1, cwf_hnd = -1, i, j, k, percent = 0, old_percent = -1, kill_count = 0, prev_point, next_point, val_cnt,
                     strongest = 0, last_valid = 0, next_valid = 0;
  float              strongest_amp, amp, next_amp, adiff, tdiff, lbounds = 20.0, ubounds = 40.0;
  CZMIL_CPF_Header   cpf_header;
  CZMIL_CPF_Data     cpf;
  CZMIL_CWF_Header   cwf_header;
  CZMIL_CWF_Data     cwf;
  uint8_t            channel[9] = {NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse, NVFalse}, check_it, mod_flag;
  char               c;
  extern char        *optarg;
  extern int         optind;


  fprintf (stderr, "\n\n %s \n\n\n", VERSION);


  while ((c = getopt (argc, argv, "d:L:U:12345679")) != EOF)
    {
      switch (c)
        {
        case 'd':
          sscanf (optarg, "%f", &adiff);
          break;

        case 'L':
          sscanf (optarg, "%f", &lbounds);
          break;

        case 'U':
          sscanf (optarg, "%f", &ubounds);
          break;

        case '1':
          channel[0] = NVTrue;
          break;

        case '2':
          channel[1] = NVTrue;
          break;

        case '3':
          channel[2] = NVTrue;
          break;

        case '4':
          channel[3] = NVTrue;
          break;

        case '5':
          channel[4] = NVTrue;
          break;

        case '6':
          channel[5] = NVTrue;
          break;

        case '7':
          channel[6] = NVTrue;
          break;

        case '9':
          channel[8] = NVTrue;
          break;

        default:
          usage ();
          break;
        }
    }


  /*  Make sure we got the mandatory file name argument and at least one channel to filter.  */

  check_it = NVFalse;

  for (i = 0 ; i < 9 ; i++)
    {
      if (channel[i])
        {
          check_it = NVTrue;
          break;
        }
    }


  if (!check_it || optind >= argc) usage ();


  strcpy (cpf_file, argv[optind]);


  if (!strstr (cpf_file, ".cpf")) usage ();


  if ((cpf_hnd = czmil_open_cpf_file (cpf_file, &cpf_header, CZMIL_UPDATE)) < 0)
    {
      czmil_perror ();
      exit (-1);
    }


  strcpy (cwf_file, cpf_file);
  sprintf (&cwf_file[strlen (cwf_file) - 4], ".cwf");

  if ((cwf_hnd = czmil_open_cwf_file (cwf_file, &cwf_header, CZMIL_READONLY)) < 0)
    {
      czmil_perror ();
      exit (-1);
    }


  fprintf (stderr, "\n\n File : %s\n\n", cpf_file);


  for (i = 0 ; i < cpf_header.number_of_records ; i++)
    {
      if (czmil_read_cpf_record (cpf_hnd, i, &cpf) != CZMIL_SUCCESS)
        {
          czmil_perror ();
          exit (-1);
        }


      if (czmil_read_cwf_record (cwf_hnd, i, &cwf) != CZMIL_SUCCESS)
        {
          czmil_perror ();
          exit (-1);
        }


      mod_flag = NVFalse;

      for (j = 0 ; j < 8 ; j++)
        {
          if (channel[j] && cpf.optech_classification[j] < CZMIL_OPTECH_CLASS_WATER)
            {
              check_it = NVFalse;


              /*  First, make sure there are at least 2 valid returns and save the index of the last valid return.  */

              val_cnt = 0;
              for (k = 0 ; k < cpf.returns[j] ; k++)
                {
                  if (!(cpf.channel[j][k].status & CZMIL_RETURN_INVAL))
                    {
                      val_cnt++;
                      last_valid = k;
                    }
                }


              /*  No point in checking if there are less than 2 valid returns.  */

              if (val_cnt >= 2)
                {
                  /*  Find the valid return with the largest amplitude.  */

                  strongest_amp = -9999.0;
                  for (k = 0 ; k < cpf.returns[j] ; k++)
                    {
                      if (!(cpf.channel[j][k].status & CZMIL_RETURN_INVAL))
                        {
                          prev_point = (int32_t) cpf.channel[j][k].interest_point;
                          next_point = prev_point + 1;


                          /*  Interpolate the interest point amplitude.  */

                          amp = (float) cwf.channel[j][prev_point] + ((float) cwf.channel[j][next_point] - (float) cwf.channel[j][prev_point]) * 
                            ((cpf.channel[j][k].interest_point - (float) prev_point) / ((float) next_point - (float) prev_point));


                          if (amp > strongest_amp)
                            {
                              strongest_amp = amp;
                              strongest = k;
                            }
                        }
                    }


                  /*  If the strongest return is the last valid return then we don't need to check it.  */

                  if (strongest != last_valid)
                    {
                      /*  Get the index of the next valid return.  */

                      for (k = strongest + 1 ; k < cpf.returns[j] ; k++)
                        {
                          if (!(cpf.channel[j][k].status & CZMIL_RETURN_INVAL))
                            {
                              next_valid = k;
                              break;
                            }
                        }


                      /*  If the next valid return is not the last valid return then we're gonna bail out so that we don't kill points in the canopy.  */

                      if (next_valid == last_valid)
                        {
                          /*  Check to see that the time difference between the strongest return and the next valid return is within our bounds.  */

                          tdiff = cpf.channel[j][next_valid].interest_point - cpf.channel[j][strongest].interest_point;

                          if (tdiff >= lbounds && tdiff <= ubounds)
                            {
                              /*  Get the amplitude of the after-pulse interest point.  */

                              prev_point = (int32_t) cpf.channel[j][next_valid].interest_point;
                              next_point = prev_point + 1;

                              next_amp = (float) cwf.channel[j][prev_point] + ((float) cwf.channel[j][next_point] - (float) cwf.channel[j][prev_point]) * 
                                ((cpf.channel[j][next_valid].interest_point - (float) prev_point) / ((float) next_point - (float) prev_point));


                              /*  Check to see if the difference between the strongest interest point amplitude and the after-pulse interest point
                                  amplitude meets or exceeds our difference threshold.  */

                              if (strongest_amp - next_amp >= adiff)
                                {
                                  cpf.channel[j][next_valid].status |= CZMIL_RETURN_FILTER_INVAL;
                                  cpf.channel[j][next_valid].filter_reason = CZMIL_AFTER_PULSE;
                                  kill_count++;
                                  mod_flag = NVTrue;
                                }
                            }
                        }
                    }
                }
            }
        }


      if (mod_flag)
        {
          if (czmil_update_cpf_return_status (cpf_hnd, i, &cpf) != CZMIL_SUCCESS)
            {
              czmil_perror ();
              exit (-1);
            }
        }


      percent = NINT (((float) i / (float) cpf_header.number_of_records) * 100.0);
      if (old_percent != percent)
        {
          fprintf (stdout, "%3d%% processed    \r", percent);
          fflush (stdout);
          old_percent = percent;
        }
    }


  fprintf (stdout, "100%% processed, %d invalidated\n", kill_count);
  fflush (stdout);

  czmil_close_cwf_file (cwf_hnd);
  czmil_close_cpf_file (cpf_hnd);

  return (0);
}
