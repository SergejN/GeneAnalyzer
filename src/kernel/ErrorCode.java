/*
    File:
        ErrorCode.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        GeneAnalyzer error codes.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel;


public enum ErrorCode 
{        
    Ok,
    ObjectNotFound,
    ObjectExists,
    FileDoesNotExist,
    FileExists,
    DirectoryDoesNotExist,
    DirectoryExists,
    IOError,
    CancelledByUser,
    ExecutionError,
    NoDatasetLoaded,
    SelectionIsEmpty,
    InvalidParameter,
    VMError
}
