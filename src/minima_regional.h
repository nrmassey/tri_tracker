/******************************************************************************
** Program : minima_regional.h
** Author  : Neil Massey
** Date    : 20/04/15
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded and then processed in some way.  This is a
**           virtual class and needs to be inherited from.  It was written
**           as many feature identification routines have the same functions
**           after the data processing.  So all you have to do is write the
**           data processing function (and overload anything as necessary)
******************************************************************************/