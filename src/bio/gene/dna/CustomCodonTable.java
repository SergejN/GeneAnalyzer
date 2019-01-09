/*
    File:
        CustomCodonTable.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Encapsulates a custom codon table, which can be loaded from an
        XML file.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */
package bio.gene.dna;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import kernel.ErrorCode;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;

/*
       XML File Format:

   <codontable>
       <name>Test</name>
       <codon>
          <sequence>ATG</sequence>
          <terminal>false</terminal>
          <start>true</start>
          <tlc>Met</tlc>
          <olc>M</olc>
       </codon>
       <codon>
          <sequence>CTC</sequence>
          <terminal>false</terminal>
          <start>false</start>
          <tlc>Leu</tlc>
          <olc>L</olc>
       </codon>
    </codontable>
 */
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

public class CustomCodonTable implements ICodonTable
{
    public static final String CODONTABLE_EXTENTION = ".xml";
    
    public static class Codon
    {
        private String  strAA;
        private String  strTLC;
        private String  strOLC;
        private String  strCodon;
        private int     iFold;
        private boolean bTerminal;
        private boolean bStart;

        
        public Codon(String strAA, String strTLC, String strOLC, String strCodon, boolean bTerminal, boolean bStart)
        {
            // Amino acid name.
            char[] tmp = strAA.toLowerCase().toCharArray();
            tmp[0] = Character.toUpperCase(tmp[0]);
            this.strAA = new String(tmp);
            // Three-letter-code.
            tmp = strTLC.toLowerCase().toCharArray();
            tmp[0] = Character.toUpperCase(tmp[0]);
            this.strTLC = new String(tmp);
            // One-letter-code.
            this.strOLC = strOLC.toUpperCase();
            // Codon sequence.
            this.strCodon = strCodon.substring(0,3).toUpperCase();
            iFold = 1;
            this.bTerminal = bTerminal;
            this.bStart = bStart;
        }

        /**
         *  Returns true, only if the codon is valid, i.e. consists only of three A, C, G or T.
         *
         *  @return
         */
        public boolean isValid()
        {
            return strCodon.matches("[ACGT][ACGT][ACGT]");
        }
    };
    
    
    private String strName = null;
    private HashMap<String, Codon> codons = null;
    
    /**
     *  Use a private constructor to avoid creating uninitialized 
     *  instances of CustomCodonTable.
     */
    private CustomCodonTable()
    {
        codons = new HashMap<String, Codon>();
    }
    
    public String getName()
    {
        return strName;
    }

    public boolean isTerminal(String strCodon)
    {
        Codon c = codons.get(strCodon);
        return (c!=null) ? c.bTerminal : false;
    }

    public String getAminoAcid(String strCodon, TYPE type)
    {
        Codon aa = codons.get(strCodon);
        if(aa==null)
            return null;
        else
        {
            switch(type)
            {
                case FullName:
                    return aa.strAA;
                case ThreeLetterCode:
                    return aa.strTLC;
                case OneLetterCode:
                    return aa.strOLC;
                default:
                    return null;
            }
        }
    }

    public boolean areSynonymous(String strCodon1, String strCodon2)
    {
        Codon aa1 = codons.get(strCodon1.toUpperCase());
        Codon aa2 = codons.get(strCodon2.toUpperCase());
        return (aa1!=null && aa2!=null) ? aa1.strOLC.equals(aa2.strOLC) : false;
    }

    public int getFoldFamily(String strCodon)
    {
        Codon aa = codons.get(strCodon);
        return (aa!=null) ? aa.iFold : 0;
    }    
        
    public boolean isStartCodon(String strCodon)
    {
        return false;
    }
    
    /**
     *  Creates the codons table from XML file. If the file does not exist
     *  the method returns null.
     * 
     *  If the file is an invalid codon table file or does not have 64 unique codon
     *  entries the method returns null.
     * 
     *  @param strFilename
     *  @return
     */
    public static CustomCodonTable loadFromFile(String strFilename)
    {
        if(strFilename==null || strFilename.isEmpty())
            return null;
        File f = new File(strFilename);
        if(!f.exists())
            return null;
        // Load the table.
        Vector<Codon> codons = new Vector<Codon>();
        String strName = null;
        SAXBuilder builder = new SAXBuilder();
        try
        {
            Document document = builder.build(f);
            Element rootElement = document.getRootElement();
            Element name = rootElement.getChild("name");
            strName = name.getText();
            List<Element> entries = rootElement.getChildren("codon");
            for(Element el:entries)
            {
                Element seq = el.getChild("sequence");
                Element term = el.getChild("terminal");
                Element start = el.getChild("start");
                Element aa = el.getChild("aminoacid");
                Element tlc = el.getChild("tlc");
                Element olc = el.getChild("olc");
                Codon c = new Codon(aa.getText(), tlc.getText(), olc.getText(), seq.getText(), 
                        new Boolean(term.getText()), new Boolean(start.getText()));
                if(c.isValid())
                    codons.add(c);
            }            
        }
        catch(Exception e)
        {
            return null;
        }
        // Check, if all possible codons are covered.
        if(codons.size()!=64)
            return null;
        return create(codons.toArray(new Codon[64]), strName);
    }
    
    /**
     *  Creates a new custom codon table using the specified codons. If the
     *  number of unique codons is not 64, the method returns null. If the strName
     *  is null or empty, the new table has the name "Unnamed".
     * 
     *  @param aas
     *  @param strName
     *  @return
     */
    public static CustomCodonTable create(Codon[] aas, String strName)
    {     
        if(strName==null || strName.isEmpty())
            strName = "Unnamed";
        if(aas.length!=64)
            return null;
        CustomCodonTable cct = new CustomCodonTable();
        cct.strName = strName;
        for(Codon c:aas)
            cct.codons.put(c.strCodon, c);
        if(cct.codons.size()!=64)
            return null;  
        // Calculate the fold for each codon.
        char[] bases = {'A', 'C', 'G', 'T'};
        for(Codon c:aas)
        {
            char[] tmp = c.strCodon.toCharArray();
            int nMatches = 1;
            for(char base:bases)
            {
                if(base==c.strCodon.charAt(2))
                    continue;
                tmp[2] = base;
                String s = cct.getAminoAcid(new String(tmp), TYPE.OneLetterCode);
                if(s.equalsIgnoreCase(c.strOLC))
                    nMatches++;
            }
            c.iFold = nMatches;
        }
        return cct;
    }
    
    /**
     *  Saves the codon table into the file.
     * 
     *  @param strFilename
     *  @return
     */
    public ErrorCode saveToFile(String strFilename)
    {
        Element root = new Element("codontable");
        Element name = new Element("name");
        name.setText(strName);
        root.addContent(name);
        Iterator<String> it = codons.keySet().iterator();
        while(it.hasNext())
        {        
            String strSeq = (String)it.next();
            Codon c = codons.get(strSeq);
            // Create XML entry.
            Element codon = new Element("codon");
            // Codon sequence.
            Element seq   = new Element("sequence");
            seq.setText(strSeq);
            codon.addContent(seq);
            // Start codon.
            Element start  = new Element("start");
            start.setText(Boolean.toString(c.bStart));
            codon.addContent(start);
            // Terminal codon.
            Element term  = new Element("terminal");
            term.setText(Boolean.toString(c.bTerminal));
            codon.addContent(term);
            // Amino acid.
            Element aa  = new Element("aminoacid");
            aa.setText(c.strAA);
            codon.addContent(aa);
            // TLC.
            Element tlc  = new Element("tlc");
            tlc.setText(c.strTLC);
            codon.addContent(tlc);
            // OLC.
            Element olc  = new Element("olc");
            olc.setText(c.strOLC);
            codon.addContent(olc);
            root.addContent(codon);
        }
        // Save XML.
        try
        {
            FileWriter fw = new FileWriter(new File(strFilename));
            Format format = Format.getPrettyFormat();
            new XMLOutputter(format).output(new Document(root), fw);
            fw.close();
        }
        catch(Exception e)
        {
            return ErrorCode.IOError;
        }
        return ErrorCode.Ok;
    }    
}
