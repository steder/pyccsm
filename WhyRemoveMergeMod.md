# Moving merge functions into the coupler #

merge\_mod.F90 was part of the original CPL6 coupler.  However, merge\_mod.F90 was very tightly "coupled" with the logic inside of the Coupler.  Essentially, the 2 functions in merge\_mod.F90 belonged with the coupler but were simply split out into their own file to make the coupler main.F90 file shorter and a bit easier to read.

I felt that this merge code was so specific to this specific coupler class that it needed to be included in the coupler.py file itself.

Additionally, the newest main.py scripts are extremely short and easy to read.  Considerably easier to digest then the old main.F90.
