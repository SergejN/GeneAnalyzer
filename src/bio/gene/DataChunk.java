/*
    File:
        DataChunk.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Represents a data chunk of any kind.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene;

import java.util.HashMap;


public abstract class DataChunk 
{
    /**
     *  Properties object can be used to store any user-defined data.
     */
    protected HashMap<String, Object> properties = null;
    
    
    /**
     *  Adds the property to the properties list. If the name
     *  of the property is null or empty it is not added.
     * 
     *  @param strName      property name
     *  @param value        property value
     *  @return the previous value associated with strName or null if
     *          there were no such key or new key is null or empty.
     */
    public Object addProperty(String strName, Object value)
    {
        if(strName==null || strName.isEmpty())
            return null;
        return properties.put(strName, value);
    }
    
    
    /**
     *  Returns the specified property or null if such property
     *  does not exist. If strName is null or empty, null is returned.
     * 
     *  @param strName      property name
     *  @return
     */
    public Object getProperty(String strName)
    {
        if(strName==null || strName.isEmpty())
            return null;
        return properties.get(strName);
    }
}
