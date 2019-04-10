# gestures, default fMRI scenario.

# parameters

## screen
$screen_width = 1024;
$screen_height = 768;
$bitmap_width = 900;

## scanning
$mri = 1; # fMRI mode. false = 0, true = 1.
$mri_test_by_emulation = 1; # for mri mode testing with the internal mri pulse generator.
$scan_period = 850.0;

## port output
$port_output = 0; # false = 0, true = 1.

## scheme
$picture_duration = 3000; # ms
$fixation_port_code = 2;

$scheme_files_directory = "scheme_gestures_fMRI";
### The scheme-files-directory must be located in the 'schemes' directory.
### See the README file for details.

## the base scenario
TEMPLATE "..\\..\\experiment-lib\\gestures.tem";