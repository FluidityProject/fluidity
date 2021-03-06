#include <h5partPluginInfo.h>
#include <avth5partFileFormat.h>
#include <avtMTMDFileFormatInterface.h>
#include <avtGenericDatabase.h>

// ****************************************************************************
//  Method:  h5partCommonPluginInfo::GetDatabaseType
//
//  Purpose:
//    Returns the type of a h5part database.
//
//  Programmer:  cristina -- generated by xml2info
//  Creation:    Mon Feb 27 13:53:31 PST 2006
//
// ****************************************************************************
DatabaseType
h5partCommonPluginInfo::GetDatabaseType()
{
    return DB_TYPE_MTMD;
}

// ****************************************************************************
//  Method:  h5partCommonPluginInfo::GetDefaultExtensions
//
//  Purpose:
//    Returns the default extensions for a h5part database.
//
//  Programmer:  cristina -- generated by xml2info
//  Creation:    Mon Feb 27 13:53:31 PST 2006
//
// ****************************************************************************
std::vector<std::string>
h5partCommonPluginInfo::GetDefaultExtensions()
{
    std::vector<std::string> defaultExtensions;
    defaultExtensions.push_back("h5part");

    return defaultExtensions;
}

// ****************************************************************************
//  Method: h5partCommonPluginInfo::Setuh5partDatabase
//
//  Purpose:
//      Sets up a h5part database.
//
//  Arguments:
//      list    A list of file names.
//      nList   The number of timesteps in list.
//      nBlocks The number of blocks in the list.
//
//  Returns:    A h5part database from list.
//
//  Programmer: cristina -- generated by xml2info
//  Creation:   Mon Feb 27 13:53:31 PST 2006
//
// ****************************************************************************
avtDatabase *
h5partCommonPluginInfo::SetupDatabase(const char *const *list,
                                   int nList, int nBlock)
{
    return new avtGenericDatabase(
               new avtMTMDFileFormatInterface(
                   new avth5partFileFormat(list[0])));
}
