# gestures, default practice scenario.

# parameters

## screen
$screen_width = 1920;
$screen_height = 1080;
$bitmap_width = 536;

## scanning
$mri = 0; # fMRI mode. false = 0, true = 1.
$mri_test_by_emulation = 0; # for mri mode testing with the internal mri pulse generator.
$scan_period = 1000.0;

## port output
$port_output = 0; # false = 0, true = 1.

## scheme
$picture_duration = 3000; # ms
$fixation_port_code = 2;

$scheme_files_directory = "scheme_gestures_practice";
### The scheme-files-directory must be located in the 'schemes' directory.
### See the README file for details.

## the base scenario
TEMPLATE "..\\..\\experiment-lib\\gestures.tem";