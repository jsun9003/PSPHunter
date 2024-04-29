# Standalone tools

These standalone tools written by Jun Sun _et al_., can be used for the prediction of phase seperating properties. if
you run the scripts after activating your `pixi` environment (`pixi shell`), everying should work.

## System calls

The perl scripts use system calls to manage the pipelines.  `predict_DrivingRegion.pl` uses a call to `bedtools`, `rm -rf` and `Intest-apply-test.py`.  The other perl scripts just make use of linux system commands or the aformentioned python script.
