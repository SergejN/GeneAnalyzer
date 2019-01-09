/*
    File:
        OSRoutines.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Encapsulates a few helpful routines, which make working with the OS more comfortable.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.util.Vector;


public class OSRoutines
{
    /**
     *  Returns the array of files in the directory using the specified filename
     *  filter. The method also enters into the subdirectories.
     *
     *  Remarks:
     *      If filter is null, the files are not filtered.
     *      If dir is null or points to a non-existing directory, the method
     *      returns null.
     *      If the directory is empty or there are no files which match the filter
     *      the method returns an array of length 0.
     *
     *  @param dir
     *  @param filter
     *  @return
     */
    public static File[] listFiles(File dir, FilenameFilter filter)
    {
        if(dir==null || !dir.exists())
            return null;
        Vector<File> files = listFileNames(dir, filter);
        if(files.size()==0)
            return new File[]{};
        return files.toArray(new File[1]);
    }

    private static Vector<File> listFileNames(File dir, FilenameFilter filter)
    {
        Vector<File> files = new Vector<File>();
        File[] children = dir.listFiles(filter);
        for(File f:children)
        {
            if(f.isFile())
                files.add(f);
            else if(f.isDirectory())
            {
                Vector<File> list = listFileNames(f, filter);
                files.addAll(list);
            }
        }
        return files;
    }

    /**
     *  Copies the file.
     *  The method returns false if one of the following is true:
     *      - src or dest is null
     *      - src points to a non-existing file
     *      - copying fails due to OS error
     *
     *  @param src
     *  @param dest
     *  @return
     */
    public static boolean copyFile(File src, File dest)
    {
        if(src==null || !src.exists() || dest==null)
            return false;
        try
        {
            FileInputStream fis  = new FileInputStream(src);
            FileOutputStream fos = new FileOutputStream(dest);
            byte[] buf = new byte[1024];
            int i = 0;
            while ((i = fis.read(buf)) != -1)
                fos.write(buf, 0, i);
            fis.close();
            fos.close();
            return true;
        }
        catch (Exception e)
        {
            return false;
        }
    }
}
