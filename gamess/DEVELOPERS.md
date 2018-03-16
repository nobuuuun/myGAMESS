# For Developers
## GAMESS Fortran 90 Coding Policy

The following rules must be adhered to in making changes in GAMESS.  These rules are important in maintaining portability and stability. Code that does not follow these rules will not be permitted in GAMESS. There will be no exceptions.

Golden Rule: make sure your code not only has no compiler errors (try as many compilers as possible), but that it also has no FTNCHEK diagnostics.  FTNCHEK and an updated version for Fortran 90 by Robert Moniot are a fantastic static analyzer, and only code that passes will be accepted into GAMESS. You will make everyone’s lives easier by performing these checks BEFORE submittal to GAMESS

Silver Rule: All code submissions and revisions tracking will be done using the GAMESS repository on GitHub.  Information on how to submit new features can be found here: https://github.com/gms-bbg/gamess/wiki

Rule 1. If there is a way to write code that works on all computers, write it that way.   Putting in preprocessor statements for different types of computers should be your last resort.  If it is necessary to add lines that are specific to one computer, put in code for all other supported machines.  If you don't have access to all of the types of supported hardware, you can look at the other machine specific examples found in GAMESS, or ask for help from someone who does understand the various architectures.  If a module does not already contain some machine specific statements, be especially reluctant to introduce dependencies.

Rule 2. Write a double precision program, and let the preprocessor handle any conversion to single precision, when that is necessary:

  * Use `Implicit none`. Using `Implicit Double Precision (A-H, O-Z)` specification statements is allowed (temporarily) only for existing code, but is not allowed for new code.  `Real*8` should not be used and `Integer` is just `Integer`.
  * Floating point constants given using scientific notation should be entered as a number with a decimal point, a number after the decimal, a "D", and a 2 digit exponent.  Example `1.23D+01`
  * Double precision BLAS names are used throughout, for example `DDOT` instead of `SDOT`, and `DGEMM` instead of `SGEMM`.

Rule 3. Fortran allows for generic functions.  Thus the routine `SQRT` should be used in place of `DSQRT`, as this will automatically be given the correct precision by the compilers.  Use `ABS`, `COS`, `INT`, etc.  Your compiler manual will tell you all the generic names.

Rule 4. Every routine in GAMESS begins with a line containing the name of the module and the routine.  An example is `!*MODULE xxxxxx *DECK yyyyyyy`.  Here `xxxxxx` is the name of the module and `yyyyyyy` is the name of the routine.  This rule makes it easier to read and maneuver around the GAMESS code for everyone.

Rule 5. Whenever a change is made to a module, this should be recorded in the Git change log included with the submitted code.  The change should not be recorded at the top of the file as was previously done.

Rule 6. Do not use tab to indent.  Instead, each new layer of embedding should be indented with 4 spaces.  Free form formatting should be used but lines should extend no further than the 80th column to maintain readability.

Rule 7. Stick to the Fortran naming convention for integer (I-N) and floating point variables (A-H,O-Z) whenever possible to maintain readability.

Rule 8. Always ensure that allocated arrays are deallocated properly and avoid allocations within a loop whenever possible.

Rule 9. `GOTO` statements are to be replaced with one of the many newer structures available in Fortran 90, such as `DO` loops or `CASE` statements.  Labels for `DO` loops can be included, but are not required.

Rule 10. All open, rewind, and close operations on sequential files must be performed with the subroutines `SEQOPN`, `SEQREW`, and `SEQCLO` respectively.  You can find these routines in `IOLIB` - they are easy to use.  `SQREAD`, `SQWRIT`, and various integral I/O routines like `PREAD` are used to process the contents of such files.  The variable `DSKWRK` tells if you are processing a distributed file (one split between all compute processes, `DSKWRK=.TRUE.`) or a single file on the master process (`DSKWRK=.FALSE.`, resulting in broadcasts of the data from the master to all other CPUs).

Rule 11. All READ and WRITE statements for the formatted files 5, 6, 7 (variables `IR`, `IW`, `IP`, or named files `INPUT`, `OUTPUT`, `PUNCH`) must be performed only by the master task.  Therefore, these statements must be enclosed in `IF (MASWRK) THEN` clauses.  The `MASWRK` variable is found in the PAR module, and is true on the master process only.  This avoids duplicate output from the other processes.

Rule 12. All error termination is done by `CALL ABRT` rather than a `STOP` statement.  Since this subroutine never returns, it is ok to follow it with a `STOP` statement, as compilers may not be happy without a `STOP` as the final executable statement in a routine.  The purpose of calling `ABRT` is to make sure that all parallel tasks get shut down properly.

Rule 13. Naming: All new subroutines must start with a prefix indicating the part of the code that they are servicing. The prefix should be in the form of prefixname_subroutinename.  For example, a subroutine for calculating energy used by a Monte Carlo run would be named `MC_energy` while a subroutine for shift FMO atoms would be named `FMO_shift`.

Rule 14. Any new subroutines must use the `intent` attribute for all arguments that are passed to the subroutine. Any arguments that are not changed in the subroutine should be intent in, while any arguments that are only output should have intent out. All others can have intent inout.  For example, a subroutine `x_sub` with one argument to be changed `real A` and one that should not be changed `integer B` would include the following lines:

```
subroutine(A, B)
integer, intent(in) :: B
real, intent(out) :: A
```

# GAMESS Doxygen Rules, Guidelines, and Examples 
### (Updated: Feb. 6, 2013)
Basic use of Doxygen Tags
The Doxygen comment block goes under the `MODULE`/`DECK` line in GAMESS, but before the `SUBROUTINE` line.

All tags should begin with the comment sign `C` or `!` followed by a greater-than sign, `>`. This document uses the Fortran 90 convention of `!`. All tags begin with `@`, so the line would look like:

```
!*MODULE … DECK … 
!> @brief …
!>
!> @tagname ...
SUBROUTINE …()
```

There will be an example code at the end illustrating the proper use of the tags.

Required Tag for ALL code
This tag must be included for all modified code, even a small change in an existing subroutine including bug fixes.

  *	`@brief`
Required Tags for all code except bug fixes - recommended for bug fixes
The following tags should be used if you have additional information about the subroutine being modified. If you do not know anything other than a brief description of the subroutine, the `@brief` tag is all that is needed.

  *	`@author`
  *	`@details`
  *	`@param`, if you added a new argument to a subroutine, otherwise this is not required - it is, however, highly recommended

Required Tags for all new code
If a subroutine was not in GAMESS until this submission, the following must be included in addition to the above rules.

  *	`@param`
Optional Tags
This is a short list of some commonly used additional tags that are not required, but encouraged when applicable. You are welcome to use tags not included on this list.

  * `@note`
  * `@warning`
  * `@todo`
  * `@bug`
  * `@see`
  * `@date`
  * `@file`

Descriptions of the Tags
Here are the descriptions of the tags along with examples of use. All tags should begin with the comment sign `C` or `!` followed by a greater-than sign, `>`. This document uses the Fortran 90 convention of `!`. All tags begin with `@`, so the line would look like:

```
!> @tagname
```

The term “user” refers to anyone writing code that uses your code, not the end user of GAMESS.
`@brief`
	Provides a brief description of the subroutine. The text following this tag should only be one sentence long. Any further information should go in the details section.  This tag must be followed by a line that contains only `!>` or `C>`.
`@author`
	Take credit for your work! If you made significant contributions to a subroutine, put your name along with the date and a brief description of the contribution made. Bug fixes should use the “date” tag discussed below. If you know the original authors’ names, put them at the top of the “author” tag. If you don’t know the original author’s name, ask Mark or Mike if they know. The format should be:

```
!> @author Original Author
!> @author Your Name
!> - date- comment about change
```

For example, say a subroutine was written by John Doe in July of 2010, and edited by Jane Doe in September of 2012. Then, John Doe made an additional change in October of 2012. This would be the example comment block:

```
!> @author John Doe
!> - July, 2010- Subroutine written
!> - October, 2012- Added additional loop to account for a larger variety of cases.
!> @author Jane Doe
!> - September, 2012- Cleaned up excess code to make the routine more general.
```

`@date`
	This should be used to account for bug fixes or for changes not significant enough to warrant an “author” addition. The dates given should be in order, with the earliest being at the top. The format for the “date” tag is as follows:

```
!> @date Month, Year Author’s Name
!> - Changes made 
```

For example, if John Doe made a small bug fix with the indexing variable in October of 2012, the comment block would look like this:

```
!> @date October, 2012
!> - Bug fix with indexing variable
```

`@details`
	For additional information the “details” tag should be used. This tag should be used for any information that would be useful for a user of the subroutine. It can be as long as needed and can contain relevant information such as paper references, discussion about the algorithm, background information, etc. If following other tags, especially `@brief`, an additional line should be used to separate `@details` from `@brief` like so:

```
!> @brief This subroutine changes the coordinates from Z-Matrix to Cartesian.
!>
!> @details Since GAMESS uses Cartesian coordinates for most of the calculations, there is a 
!> need for the conversion from Z-Matrix to Cartesian coordinates. This routine is based on 
!> Model-Builder written by Mark Gordon. 
```

`@param`
	If you are editing a subroutine that was written by someone else, the `@param` tags are not required, but strongly recommended if you know what the arguments mean. Even if you only know the meaning of a sub-set of the arguments, include them with an explanation. If you add a new argument to a subroutine, you must include this tag. The format for this tag is

```
!> @param argument_1 Description of argument_1
!> @param argument_2 Description of argument_2
…
```

For example, if you had a subroutine with 4 arguments:

```
SUBROUTINE DISPLACE(COORDS,X_DISPLACE,Y_DISPLACE,Z_DISPLACE)
```
The resulting comment block would be:
```
!> @param COORDS Input coordinates
!> @param X_DISPLACE Displacement in the X direction
!> @param Y_DISPLACE Displacement in the Y direction
!> @param Z_DISPLACE Displacement in the Z direction
```

`@note`
	The `@note` tag should be used to highlight information for the user. This could be anything you would want the user to be aware of. For example, you could remind users of a particular way you intended the code to be used. The syntax is:

```
!> @note Your note
```

For example, if you made a routine that should only be called at the beginning when the input file is being read, this would be the comment block:

```
!> @note Should only be called when the input is being read in.
```

`@warning`
	Should be used when what you have to say is of more importance than a simple `@note` statement. A good use of the `@warning` tag would be: If the contents of the tag are not followed, a bug could occur. The syntax is the same as `@note`.
`@todo`
	Whenever there is additional work to be done on a subroutine, the `@todo` tag should be used. The `@todo` tag will place the line into a special list that includes all of the todo statements in one place. This tag is helpful when you need to release a subroutine with basic functionality, but you have additional functionality in mind. For example:

```
!> @todo Make this routine run in parallel.
```

`@bug`
	Any known, but not easily fixable bugs should be documented with this tag. However, this tag should not be used as a substitute for debugging or error reporting. This tag should let users know that under certain conditions a subtle bug may appear, so they can keep that in mind when they use this routine. Just like `@todo` this line is placed in a special list.

`@see`
	This tag places a “see also” section in the documentation page. It is useful for linking similar parts of the code, so users won’t have to search for it themselves. Using this tag saves people time and can help you by preventing code duplication.

`@file`
	This tag is special in that it is placed at the top of the file before any subroutine specific tags.  This is for documenting information relevant to two or more subroutines, so that the documentation isn't duplicated.  It is important that an empty line between the end of the `@file` tag and the first tag used for a subroutine.

```
!> @file information about the over all method used in the code
!>           that may involve subroutines that use very similar
!>           information and methods and so would have duplicate information.

!> @brief summary of first subroutine's function
```

Example code
Here is an example subroutine that uses many of the tags previously discussed. This will give you an idea of how Doxygen can be incorporated into GAMESS.

```
!*MODULE MATH DECK MULTIPLIER
!> @brief Multiplies two arrays and returns the result
!>
!> @detail Takes in 2 3x3 arrays and multiplies
!> by multiplying the rows of Array1 with 
!> the columns of Array2
!>
!> @author Random Grad
!> - September, 1985
!> @date January, 2012- Another Grad
!> - Bug Fix
!> @date October, 2012- John Doe
!> - Different bug fix
!> @param Array1 First array
!> @param Array2 Second array
!> @param Out Returned array
!> @todo write the rest of the subroutine
Subroutine Multiplier(Array1,Array2,Out)
    Implicit none
    Double Precision, Intent(in) :: Array1(3,3),Array2(3,3)
    Double Precision, Intent(out) :: Out(3,3)
```

# GAMESS FTN90 Requirements
Codes that contributors wish to put into GAMESS must adhere to FORTRAN standards, so that they can be compiled on commonly used computer hardware and software, before they are accepted into the git repository and eventually into a GAMESS release. This is implicit in any agreement to add codes to GAMESS. The final arbiters in this regard will be the GAMSS development team.

# Additional Resources
Also see the [GAMESS DOCUMENTATION](http://www.msg.ameslab.gov/GAMESS/documentation.html) for information on Programmer's References and Hardware Specfics.


