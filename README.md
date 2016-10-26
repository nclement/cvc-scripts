= Nathan's CVC Scripts = 
Here are some cool scripts I've found useful. I've tried to make them usable on
many different platforms by creating a Makefile.def file that houses all the
paths for various files and folders. Inside this file you'll also find some
paths that currently work on Stampede. You should just be able to change the
paths inside this file and it should work across all the scripts.

== Note to Disgruntled User == 
Please note: if you should choose to copy/move one of these scripts into a
different folder, you should know that the script will likely break. Most of
these are meant to be run from another folder. Something like:

```bash
export SCRIPT_DIR=/path/to/script/dir
$SCRIPT_DIR/getTotalEnergy.sh protein.pdb
```
