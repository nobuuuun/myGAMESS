#ifndef DDI_RUNTIME_H
#define DDI_RUNTIME_H

/**
   @file
   @brief DDI version and runtime information.
 */

/**
   @brief Maximum length of build string.
*/
#define DDI_MAX_BUILD 256

/**
   @brief Gets DDI version and subversion.
   @param[out] version DDI version.
   @param[out] subversion DDI subversion.
*/
void DDI_Get_version(int *version, int *subversion);

/**
   @brief Gets DDI build information.
   @details Gets information about DDI build, such as SMP support, communication
   implementation, and I/O support.
   @warning Build string should be at least DDI_MAX_BUILD characters long.
   @param[out] build DDI build string.
*/
void DDI_Get_build(char *build);

/** 
    @brief Prints POSIX runtime as provided by uname() function 
           and DDI version and build.
    @param[in] stream File stream.
 */
void DDI_POSIX_Runtime_print(FILE *stream);

#endif /* DDI_RUNTIME_H */
