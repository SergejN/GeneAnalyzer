/*
    File:
        ViewLoader.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        GeneAnalyzer view loader.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package plugin.loader;

import java.io.File;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.jar.Manifest;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import kernel.Main;
import plugin.classes.AGeneAnalyzerView;

/*
 *  REMARKS:
 * 
 *  Every view plugin should have two entries in the manifest file:
 *      - HostMajorVersion  minimal required host major version
 *      - HostMinorVersion  minimal required host minor version.
 *                          The loading process fails if the host tries to load
 *                          the plugin which is designed for a higher host version.
 *                          The versions are backwards compatible.
 */
public final class ViewLoader extends ClassLoader
{
    // Manifest file entries.
    private static final String HOST_MAJOR_V = "HostMajorVersion";
    private static final String HOST_MINOR_V = "HostMinorVersion";
    private static final String MAIN_CLASS   = "ViewMain";
    
    // Binary class data.
    private HashMap<String, byte[]> data;
    // Manifest.
    Manifest man;
    
    /**
     *  Constructor.
     */
    public ViewLoader()
    {
        data = new HashMap<String, byte[]>();
        man = new Manifest();
    }
    
    /**
     *  Loads the class with specified name.
     * 
     *  @param strClassName
     *  @return
     *  @throws java.lang.ClassNotFoundException
     */
    @Override
    public Class loadClass(String strClassName) throws ClassNotFoundException
    {
        return loadClass(strClassName, true);
    }

    /**
     *  Loads the class with specified name.
     * 
     *  @param strClassName
     *  @param bResolve
     *  @return
     *  @throws java.lang.ClassNotFoundException
     */
    @Override
    public Class loadClass(String strClassName, boolean bResolve) throws ClassNotFoundException
    {
        // Check the cache of the parent ClassLoader
        Class c = findLoadedClass(strClassName);
        if (c != null)
        {
            return c;
        }
        // If the class was not found in the cache try to look in system classes.
        try
        {
            c = findSystemClass(strClassName);
            if (c != null)
            {
                return c;
            }
        }
        catch (Exception e)
        {
        }
        // If class was still not found, try to look in hash map.
        if(c==null)
            c = loadClassFromHashMap(strClassName, bResolve);
        // If the class was not found in the hash map throw an exception.
        if(c==null)
            throw new ClassNotFoundException("Class "+strClassName+" could not be loaded.");
        // Otherwise return the class.
        return c;
    }

    /**
     *  Loads the class bytes from the hash map.
     * 
     *  @param strClassName
     *  @param bResolve
     *  @return
     *  @throws java.lang.ClassNotFoundException
     */
    private Class loadClassFromHashMap(String strClassName, boolean bResolve) throws ClassNotFoundException
    {
        // In class name replace all "." by "/".
        byte[] classbytes = data.get(strClassName.replaceAll("\\.", "/") + ".class");
        Class c;
        if (classbytes != null)
        {
            c = defineClass(strClassName,
                            classbytes,
                            0,
                            classbytes.length);
            if(bResolve)
                resolveClass(c);
            return c;
        }        
        else
        {
            return null;
        }
    }
    
    /**
     *  Loads the specified plugin.
     *  If one of the following is true, the method returns null:
     *      - Specified file does not exist or does not specify a file
     *      - Host-application version is not supported by the plugin
     *      - Plugin file does not have a PluginMain class
     *      - I/O error occurs
     *
     *  @param strPlugin
     *  @return
     */
    public AGeneAnalyzerView loadPlugin(String strPlugin)
    {
        // Check whether the file exists.
        File f = new File(strPlugin);
        if(!f.isFile() || !f.exists())
            return null;
        // Load JAR.
        try
        {
            ZipFile zf = new ZipFile(strPlugin);
            // Get the list of entries.
            Enumeration entries = zf.entries();
            while (entries.hasMoreElements())
            {
                ZipEntry ze = (ZipEntry) entries.nextElement();
                InputStream is = zf.getInputStream(ze);
                // If the entry is the manifest, update man member.
                if(ze.getName().endsWith("MANIFEST.MF"))
                {
                    man = new Manifest(is);
                    // Check the version.
                    int iMajor = Integer.parseInt(man.getMainAttributes().getValue(HOST_MAJOR_V));
                    if(iMajor>Main.HOST_MAJOR_VERSION)
                        return null;
                    // If the major version is the same as the major version of the host-application,
                    // check the minor version. Otherwise no check is necessary, because of the backwards
                    // compatibility.
                    else if(iMajor==Main.HOST_MAJOR_VERSION)
                    {
                        int iMinor = Integer.parseInt(man.getMainAttributes().getValue(HOST_MINOR_V));
                        if(iMinor>Main.HOST_MINOR_VERSION)
                            return null;
                    }
                }
                else
                {
                    // Otherwise read the whole JAR file entry into an array.
                    byte[] b = new byte[(int) ze.getSize()]; 
                    // Read the data from stream.
                    int nRead = 0;
                    int nBytesChunk = 0;
                    int nTotal = (int)ze.getSize();
                    while(nTotal-nRead>0)
                    {
                        nBytesChunk = is.read(b, nRead, nTotal-nRead);
                        if(nBytesChunk==-1)
                            break;
                        nRead+=nBytesChunk;
                    }                
                    data.put(ze.getName(), b);
                }
            }
            zf.close();
            // Load plugin main class.
            Class c = loadClassFromHashMap(MAIN_CLASS, true);
            if(c==null)
            {
                return null;
            }
            AGeneAnalyzerView view = (AGeneAnalyzerView)c.newInstance();
            return view;
        }
        catch (Exception e)
        {
            return null;
        }
    }
}
