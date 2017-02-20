--------------------------------------------------------------------------------

Description: A revision of Wisecondor newref.py and test.py.
             - Added detection of CNV on chromosome X and Y.
             - Improved efficiency and stability.

--------------------------------------------------------------------------------

usage: newref.py [-h] [-female] [-ignore IGNORE] [-maxbin1 MAXBIN1]
                 [-maxbin2 MAXBIN2] [-refmaxval REFMAXVAL]
                 [-refmaxrep REFMAXREP]
                 refin refout

Create a new reference table from a set of reference samples, outputs table as
pickle to a specified output file

positional arguments:
  refin                 directory containing samples (.correct) to be used as
                        reference
  refout                reference table output, used for sample testing
                        (pickle)

optional arguments:
  -h, --help            show this help message and exit
  -female               turn on if gender is female (default: False)
  -ignore IGNORE        ignore x highest scoring sample samples per bin
                        distance calculation, 0 to keep all, use if you are
                        unsure about what samples have aberrations and assume
                        at most x have the same aberration, beware: results
                        become less accurate (default: 0)
  -maxbin1 MAXBIN1      maximum number of reference bins for each target bin
                        before removing neighboring bins (default: 250)
  -maxbin2 MAXBIN2      maximum number of reference bins for each target bin
                        after removing neighboring bins (final) (default: 100)
  -refmaxval REFMAXVAL  start cutoff value for determining good quality
                        reference bins (default: 1000000)
  -refmaxrep REFMAXREP  amount of improval rounds for determining good quality
                        reference bins (default: 3)

--------------------------------------------------------------------------------

usage: test.py [-h] [-female] [-maxrounds MAXROUNDS] [-refminbin REFMINBIN]
               [-refmaxbin REFMAXBIN] [-window WINDOW]
               sample reference outfile

Calculate z-scores

positional arguments:
  sample                sample to be tested (.correct)
  reference             reference table used for within sample comparison
                        (pickle)
  outfile               output results to this file

optional arguments:
  -h, --help            show this help message and exit
  -female               turn on if gender is female (default: False)
  -maxrounds MAXROUNDS  maximum amount of rounds used to calculate z-score
                        (default: 5)
  -refminbin REFMINBIN  minimum number of reference bins, ignore target bin if
                        there are less reference bins available (default: 10)
  -refmaxbin REFMAXBIN  maximum number of reference bins, ignore any reference
                        bin after (default: 100)
  -window WINDOW        window size for sliding window approach, number of
                        bins is considered in each direction (i.e. using 3
                        results in using 3+1+3=7 bins per call) (default: 5)