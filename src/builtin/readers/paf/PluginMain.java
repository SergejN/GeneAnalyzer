/*
    File:
        PluginMain.java
 *   
    Revision:
        1.2.0.1
 * 
    Description:
        Parses the data files in Peter Andolfatto's FASTA format.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.readers.paf;

import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;
import plugin.classes.ADatasetImporter;


public class PluginMain extends ADatasetImporter
{
    private String strLastError = "";
    private int index           = 0;

    @Override
    public String GetFileExtension()
    {
        return "fa";
    }

    @Override
    public String GetFileDescription()
    {
        return "Files in Peter Andolfatto's FASTA format";
    }

    @Override
    public Dataset ImportDataset(File[] files, String strParams)
    {
        index = 0;
        Dataset ds = new Dataset();
        for(File f:files)
        {
            Dataset tmp = createDataset(f);
            ds.merge(tmp);
        }
        return ds;
    }

    public ErrorCode Initialize(AInitData initdata)
    {
        return ErrorCode.Ok;
    }

    public String GetMenuItemName()
    {
        return "Peter Andolfatto's FASTA";
    }

    public String GetName()
    {
        return "Peter Andolfatto's FASTA files reader";
    }

    public String GetDescription()
    {
        return "Reads the data in Peter Andolfatto's FASTA format";
    }

    public boolean SupportsMissingData()
    {
        return true;
    }

    public boolean SupportsAmbiguousData()
    {
        return true;
    }

    public String GetParamString()
    {
        return "";
    }
    
    public String GetLastError()
    {
        return strLastError;
    }

    
    private Dataset createDataset(File file)
    {
        index++;
        try
        {
            BufferedReader in = new BufferedReader(
                                    new InputStreamReader(
                                        new FileInputStream(file)));
            Dataset ds = new Dataset();
            GeneEntry ge = null;
            
            String strHeader = null;
            StringBuffer seq = new StringBuffer();
            
            // Contains the regions annotations. Is also used as a flag to
            // specify whether the sequence or annotation data is being parsed.
            StringBuffer annot = null;
            
            String strLine;
            while((strLine = in.readLine())!=null)
            {
                // A new header line candidate is found.
                if(strLine.startsWith(">"))
                {
                    // Create a new StrainEntry instance from the data loaded previously.
                    if(strHeader!=null)
                    {
                        GeneEntry tmp = addStrainEntry(strHeader,seq.toString(), ge);
                        if(ge!=tmp)
                        {
                            ds.addGene(ge); 
                            ge = tmp;
                        }
                    }
                    // Set new header line.
                    strHeader = strLine;
                    seq = new StringBuffer();
                }
                // Annotation part marker found. The parser does not expect
                // to find any annotation right after this mark in the same line.
                // Any SEQUENCE lines after this mark will be ignored.
                else if(strLine.startsWith(";"))
                {
                    // Make sure it is the first mark occurence. If it is not
                    // simply ignore it.
                    if(annot==null)
                    {
                        annot = new StringBuffer();
                        // Add last gene entry and last strain entry.
                        if(strHeader!=null)
                        {
                            ds.addGene(ge);
                            ds.addGene(addStrainEntry(strHeader, seq.toString(), ge));
                        }
                    }
                }
                // Annotation data.
                else if(annot!=null)
                {
                    annot.append(strLine);
                    annot.append(";");
                }
                // Sequence data.
                else
                    seq.append(strLine);
            }       
            // Offset indicates the shift caused by used ORF.
            int offset = 0;
            Pattern p = Pattern.compile("FRAME\\s+(\\d+)", Pattern.CASE_INSENSITIVE);
            Matcher m = p.matcher(annot);
            if(m.find())
                offset = Integer.parseInt(m.group(1));
            // Parse the annotation lines.
            p = Pattern.compile("([^;\\s]+)\\s+(\\d+)\\s*-\\s*(\\d+)");
            m = p.matcher(annot);
            Vector<GeneRegion> prec = new Vector<GeneRegion>(); // Precursors.
            while(m.find())
            {
                String strType = m.group(1);
                // Find the known region abbreviations.
                if(strType.equalsIgnoreCase("CDS"))
                    strType = GeneRegion.EXON;
                if(strType.equalsIgnoreCase("intron"))
                    strType = GeneRegion.INTRON;
                if(strType.equalsIgnoreCase("three_UTR"))
                    strType = GeneRegion.UTR3;
                if(strType.equalsIgnoreCase("five_UTR"))
                    strType = GeneRegion.UTR5;
                if(strType.equalsIgnoreCase("intergenic"))
                    strType = GeneRegion.INTERGENIC;
                // Create a precursor region.
                GeneRegion r = new GeneRegion(strType);
                r.setStart(Integer.parseInt(m.group(2)));
                r.setEnd(Integer.parseInt(m.group(3)));
                prec.add(r);
            }
            // Create the regions for loaded genes.
            for(int i=0;i<ds.getGenesCount();i++)
            {
                for(int j=0;j<ds.getGeneEntry(i).getStrainsCount();j++)
                {                    
                    StrainEntry tmp = ds.getGeneEntry(i).getStrainEntry(j);
                    String strSeq = tmp.getRegion(0).getSequence();
                    tmp.removeRegion(0);
                    for(int n=0;n<prec.size();n++)
                    {
                        GeneRegion pr = prec.get(n);
                        GeneRegion r = new GeneRegion(pr.getType());
                        r.setStart(pr.getStart());
                        r.setEnd(pr.getEnd());                        
                        r.setSequence(strSeq.substring(r.getStart()-1, r.getEnd()));
                        tmp.addRegion(r);
                    }
                    // If the frame is not 0, remove the frame-shift from the first exon and
                    // update start/end positions of all regions after the first exon.
                    if(offset>0)
                    {
                        boolean b = false;  // Indicates whether the first exon was found.
                        for(int n=0;n<tmp.getRegionsCount();n++)
                        {
                            GeneRegion r = tmp.getRegion(n);
                            if(b)   // First exon already found, update positions.
                            {
                                r.setStart(r.getStart()-offset);
                                r.setEnd(r.getEnd()-offset);
                            }
                            // First exon not found yet.
                            else if(tmp.getRegion(n).hasType(GeneRegion.EXON))
                            {
                                b = true;                                
                                // Only update the end position.
                                r.setEnd(r.getEnd()-offset);
                                // Remove extra bases.
                                r.setSequence(r.getSequence().substring(offset));
                            }
                        }
                    }
                }
            }
            return ds;
        }
        // Return null on error and set the error description.
        catch(FileNotFoundException e)
        {
            strLastError = String.format("File %s was not found", file.getAbsolutePath());
            return null;
        }
        catch(IOException e)
        {
            strLastError = String.format("An I/O error occured while reading", file.getAbsolutePath());
            return null;
        }
    }
    
    
    /**
     *  Parses the header line and the sequence and creates a new strain entry.
     *  After that adds this entry to the specified gene entry. If the gene entry
     *  is null, it is created first. If the parsed common gene name is different
     *  than the one of the gene entry, a new gene entry is created.
     *  The method returns the gene entry it added the strain to.
     * 
     *  @param strHeaderLine
     *  @param strSequence
     *  @param ge
     *  @return 
     */
    private GeneEntry addStrainEntry(String strHeaderLine, String strSequence, GeneEntry ge)
    {
        // Create a new strain entry.
        StrainEntry se = new StrainEntry(null, null);
        se.setChromosome("X");
        // Since the coordinates are just at the end of the file,
        // add the complete sequence as "Unnamed region" and delete it later on.
        GeneRegion reg = new GeneRegion(GeneRegion.UNNAMED);
        reg.setSequence(strSequence);
        se.addRegion(reg);
         
        // The naming pattern of PAF formatted files is as follows:
        // Every header line begins with a > symbol. This can be followed
        // by either a species name (e.g. Dsim, Dyak), or a strain name
        // (e.g. ZW190, ZH23) whereas Drosophila melanogaster is assumed to be
        // the species name or a line without any CG number ending with a "genome" - 
        // this is a genomic sequence from a DB then.
        // The current version of the parser extracts the strain and the species
        // names as well as the DB sequence. The remainder of the header line
        // is not parsed even if containing meaningful information. It is saved
        // under the name RAWHEADER in the strain's annotation.        
        // Parse the strain name.
        Pattern p = Pattern.compile(">([^_]+)", Pattern.CASE_INSENSITIVE);
        Matcher m = p.matcher(strHeaderLine);
        if(m.find())
        {
            // D.melanogaster with ZW in the strain name.
            if(m.group(1).matches("[Zz][Ww]\\d+"))
            {
                se.addPopulations(new String[]{"Africa", "Drosophila melanogaster"});
                se.setSpeciesName("Drosophila melanogaster");
                se.setStrainName(m.group(1).replace("ZW", "ZBMEL"));
            }
            // D.melanogaster with ZBMEL in the strain name.
            else if(m.group(1).matches("ZBMEL\\d+"))
            {
                se.addPopulations(new String[]{"Africa", "Drosophila melanogaster"});
                se.setSpeciesName("Drosophila melanogaster");
                se.setStrainName(m.group(1));
            }
            // D.melanogaster of any other strain.
            else if(m.group(1).matches("[Zz][HhSs]\\d+"))
            {
                se.addPopulations(new String[]{"Melanogaster (ZH/ZS)", "Drosophila melanogaster"});
                se.setSpeciesName("Drosophila melanogaster");
                se.setStrainName(m.group(1));
            }
            // D.simulans
            else if(m.group(1).startsWith("Dsim") || m.group(1).startsWith("dsim") || m.group(1).startsWith("MD99"))
            {
                se.addPopulations(new String[]{"Drosophila simulans"});
                se.setSpeciesName("Drosophila simulans");
                se.setStrainName("sim_pa");
            }
            // D.yakuba
            else if(m.group(1).startsWith("Dyak") || m.group(1).startsWith("dyak"))
            {
                se.addPopulations(new String[]{"Drosophila yakuba"});
                se.setSpeciesName("Drosophila yakuba");
                se.setStrainName("yak_pa");
            }
            // D.sechellia
            else if(m.group(1).startsWith("Dsec") || m.group(1).startsWith("dsec"))
            {
                se.addPopulations(new String[]{"Drosophila sechellia"});
                se.setSpeciesName("Drosophila sechellia");
                se.setStrainName("sec_pa");
            }
            // chrX. In this case no CG number can be parsed - add the strain
            // without comparing the gene names.
            else if(m.group(1).startsWith("chrX"))
            {
                se.addPopulations(new String[]{"Genome sequence"});
                se.setSpeciesName("Genome");
                se.setStrainName("genome_pa");
            }
            // Any other name.
            else
            {
                se.addPopulations(new String[]{m.group(1)});
                se.setSpeciesName(m.group(1));
                se.setStrainName(m.group(1));
            }
        }        
        // Parse CG number.
        p = Pattern.compile("(CG\\d+)");
        m = p.matcher(strHeaderLine);
        if(m.find())
        {
            // If the CG number is not the same as in previous entries,
            // create a new gene entry.
            if( (ge==null) || !ge.getCommonName().matches(m.group(1)+"_\\d+"))
                ge = new GeneEntry(m.group(1)+"_"+Integer.toString(index), "");            
        }
        // Add strain entry.
        ge.addStrain(se);
        // Annotate raw header.
        se.addProperty("RAWHEADER", strHeaderLine);
        return ge;
    }
}