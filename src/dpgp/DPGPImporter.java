/*
    File:
        DPGPImporter.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Imports DPGP data. The required data are *.VMA files of one chromosome
        and a *.GFF file for this chromosome.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package dpgp;

import dpgp.DatasetBuilder.Fragment;
import dpgp.gui.MainDlg;
import bio.gene.Dataset;
import gui.IWaitDialog;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;


public class DPGPImporter
{
    private String strErr = null;
    private AInitData id  = null;
    private ErrorCode lec;          // Last error code.

    /**
     *  Initializes the DPGPImporter instance.
     *
     *  @param initdata     initialization data
     *  @return
     */
    public ErrorCode initialize(AInitData initdata)
    {
        id = initdata;
        return ErrorCode.Ok;
    }

    /**
     *  Returns the backbone of the parameter string as follows:
     *      - vma       list of semicolon-separated *.VMA file names
     *      - cutoff    cutoff threshold
     *      - chr       chromosome
     *      - specs     species
     *      - strains   strains
     *      - gff       GFF file name
     *      - frags     fragments file name
     *
     *  Remarks:
     *  The last two parameters are mutually exclusive and, thus, should not be both
     *  used in the same parameter string.
     *
     *  @return
     */
    public static String getParamString()
    {
        return "vma='<FILE1>;<FILE2>' cutoff='<CUTOFF>' chr='<CHROMOSOME>' specs='<SPEC1>;<SPEC2>' "+
               "strains='<STRAIN1>;<STRAIN2>' gff='<GFFFILE>' frags='<FRAGMENTSFILE>'";
    }

    /** 
     *  Shows a dialog allowing the user to specify all parameters and then
     *  creates a dataset.
     * 
     *  @return
     */
    public Dataset importDataset()
    {
        MainDlg dlg = new MainDlg(id.sm);
        MainDlg.Options opt = dlg.getImportOptions();
        if(opt==null || (opt.gff==null && opt.frags==null && opt.fragfile==null))
        {
            lec = ErrorCode.CancelledByUser;
            return null;
        }
        if(opt.gff!=null)
            return importDataset(opt.vma, opt.gff, opt.specs, opt.strains, opt.strChromosome, opt.iCutoff);
        else if(opt.fragfile!=null)
            return importDataset(opt.vma, opt.specs, opt.strains, opt.strChromosome, opt.iCutoff, opt.fragfile);
        else
            return importDataset(opt.vma, opt.frags, opt.specs, opt.strains, opt.strChromosome, opt.iCutoff);
    }

    /**
     *  Parses the options from the string and builds the dataset.
     *
     *  Parameter string format:
     *  vma='<FILE1>;<FILE2>' cutoff='<CUTOFF>' chr='<CHROMOSOME>' specs='<SPEC1>;<SPEC2>'
     *  strains='<STRAIN1>;<STRAIN2>' gff='<GFFFILE>' frags='<FRAGMENTSFILE>'
     *
     *  @param strParams
     *  @return
     */
    public Dataset importDataset(String strParams)
    {
        if(strParams==null || strParams.isEmpty())
            return importDataset();
        else
        {
            Pattern p = Pattern.compile("vma='(.+)'\\s+"+               // 1
                                        "cutoff='(\\d+)'\\s+"+          // 2
                                        "chr='(.+)'\\s+"+               // 3
                                        "specs='(.+)'\\s+"+             // 4
                                        "strains='(.+)'\\s+"+           // 5
                                        "(gff='.+')|(frags='.+')",      // 6
                                        Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
            Matcher m = p.matcher(strParams);
            if(m.find())
            {
                File[] vma = null;
                int iCutoff = 0;
                String strChr = null;
                String[] specs = null;
                String[] strains = null;
                // VMA files.
                Vector<File> vmaFiles = new Vector<File>();
                String[] tmp = m.group(1).split("\\s*;\\s*");
                for(String s:tmp)
                {
                    File f = new File(s);
                    if(f.exists())
                        vmaFiles.add(f);
                }
                if(vmaFiles.size()==0)
                    return importDataset();
                else
                    vma = vmaFiles.toArray(new File[1]);
                // Cutoff.
                iCutoff = Integer.parseInt(m.group(2));
                if(iCutoff<1)
                    return importDataset();
                // Chromosome.
                strChr = m.group(3);
                // Species and strains.
                specs = m.group(4).split("\\s*;\\s*");
                strains = m.group(5).split("\\s*;\\s*");
                // GFF file.
                if(m.group(6).startsWith("gff"))
                {
                    Pattern p2 = Pattern.compile("='([^']+)'",
                            Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
                    Matcher m2 = p2.matcher(m.group(6));
                    if(m2.find())
                    {
                        File f = new File(m2.group(1));
                        if(f.exists())
                            return importDataset(vma, f, specs, strains, strChr, iCutoff);
                        else
                            return importDataset();
                    }
                    else
                        return importDataset();
                }
                // Fragments file.
                else if(m.group(6).startsWith("frags"))
                {
                    Pattern p2 = Pattern.compile("='([^']+)'",
                            Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE);
                    Matcher m2 = p2.matcher(m.group(6));
                    if(m2.find())
                    {
                        File f = new File(m2.group(1));
                        if(f.exists())
                            return importDataset(vma, specs, strains, strChr, iCutoff, f);
                        else
                            return importDataset();
                    }
                    else
                        return importDataset();
                }
                else
                    return importDataset();
            }
            else
                return importDataset();
        }
    }

    /**
     *  Parses VMA files and the GFF file and returns the dataset.
     *
     *  @param vma
     *  @param gff
     *  @param specs
     *  @param strains 
     *  @param strChromosome
     *  @param iCutoff
     *  @return
     */
    public Dataset importDataset(File[] vma, File gff, String[] specs, String[] strains, String strChromosome, int iCutoff)
    {
        id.wd.show(IWaitDialog.TYPE.Import);
        // First, create the multiple alignment.
        DPGPAlignment ali = buildAlignment(vma, specs, strains, iCutoff);
        if(ali==null)
        {
            id.wd.close();
            lec = ErrorCode.ExecutionError;
            return null;
        }
        // Second, parse GFF file and create fragments.
        GFFParser gffp = new GFFParser();
        HashMap<String, Vector<DatasetBuilder.Fragment>> fragments = gffp.createFragments(gff, strChromosome);
        if(fragments==null)
        {
            id.wd.close();
            strErr = gffp.getLastErrorString();
            lec = ErrorCode.ExecutionError;
            return null;
        }
        // Finally, create the dataset.
        return buildDataset(ali, fragments, strChromosome);
    }

    /**
     *  Parses VMA files and creates a dataset using the provided fragments.
     *
     *  @param vma
     *  @param fragments
     *  @param specs
     *  @param strains
     *  @param strChromosome
     *  @param iCutoff
     *  @return
     */
    public Dataset importDataset(File[] vma, HashMap<String, DatasetBuilder.Fragment> fragments, String[] specs, String[] strains, String strChromosome, int iCutoff)
    {
        id.wd.show(IWaitDialog.TYPE.Import);
        // First, create the multiple alignment.
        DPGPAlignment ali = buildAlignment(vma, specs, strains, iCutoff);
        if(ali==null)
        {
            id.wd.close();
            lec = ErrorCode.ExecutionError;
            return null;
        }
        // Then, create a separate gene for each fragment.
        HashMap<String, Vector<DatasetBuilder.Fragment>> frags = new HashMap<String, Vector<Fragment>>();
        Iterator<String> it = fragments.keySet().iterator();
        while(it.hasNext())
        {
            String strName = it.next();
            Vector<Fragment> fr = new Vector<Fragment>();
            fr.add(fragments.get(strName));
            frags.put(strName, fr);
        }
        // Create the dataset.
        return buildDataset(ali, frags, strChromosome);
    }

    /**
     *  Parses VMA files and creates a dataset using the provided fragments file.
     *
     *  Remarks:
     *      Note different parameter ordering in this method compared to the
     *      one using GFF file. This method does not parse GFF file, but instead
     *      parses a plain text file with fragments coordinates.
     *
     *  @param vma
     *  @param specs
     *  @param strains
     *  @param strChromosome
     *  @param iCutoff
     *  @param fragments
     *  @return
     */
    public Dataset importDataset(File[] vma, String[] specs, String[] strains, String strChromosome, int iCutoff, File fragments)
    {
        id.wd.show(IWaitDialog.TYPE.Import);
        // First, create the multiple alignment.
        DPGPAlignment ali = buildAlignment(vma, specs, strains, iCutoff);
        if(ali==null)
        {
            id.wd.close();
            lec = ErrorCode.ExecutionError;
            return null;
        }
        HashMap<String, Vector<DatasetBuilder.Fragment>> frags = new HashMap<String, Vector<Fragment>>();
        // Read the file.
        try
        {
            BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(fragments)));
            String strLine;
            while((strLine = in.readLine())!=null)
            {
                if(strLine.matches("[^;]+;[^;]+;\\d+;\\d+"))
                {
                    String[] tmp = strLine.split("\\s*;\\s*");
                    String strName = tmp[0];
                    String strType = tmp[1];
                    int iStart = Integer.parseInt(tmp[2]);
                    int iEnd = Integer.parseInt(tmp[3]);
                    Vector<DatasetBuilder.Fragment> vf = new Vector<Fragment>();
                    DatasetBuilder.Fragment f = new Fragment();
                    f.bLeading = (iStart<=iEnd);
                    f.strType = strType;
                    f.iStart = Math.min(iStart, iEnd);
                    f.iEnd = Math.max(iStart, iEnd);
                    vf.add(f);
                    frags.put(strName, vf);
                }
            }
        }
        catch(FileNotFoundException e)
        {
            id.wd.close();
            strErr = "File not found";
            lec = ErrorCode.FileDoesNotExist;
            return null;
        }
        catch(IOException e)
        {
            id.wd.close();
            strErr = "An I/O error occured while parsing the file";
            lec = ErrorCode.IOError;
            return null;
        }
        if(frags.size()==0)
        {
            id.wd.close();
            strErr = "No valid entries found in the fragments file";
            lec = ErrorCode.ExecutionError;
            return null;
        }
        // Create the dataset.
        return buildDataset(ali, frags, strChromosome);
    }

    /**
     *  Builds the multiple alignment using the provided VMA files. THe VMA
     *  files should represent different species but the same chromosome, otherwise
     *  the alignment is incorrect.
     *  
     *  Remarks:
     *  If one of the following is true, the method returns null:
     *      - vma is null
     *      - specs is null
     *      - strains is null
     *      - vma, specs and strains have different lengths
     *
     *  @param vma
     *  @param specs
     *  @param strains
     *  @param iCutoff
     *  @return
     */
    public DPGPAlignment buildAlignment(File[] vma, String[] specs, String[] strains, int iCutoff)
    {
        if(vma==null || specs==null || strains==null || vma.length!=specs.length || vma.length!=strains.length)
        {
            strErr = "Invalid parameter";
            lec = ErrorCode.InvalidParameter;
            return null;
        }        
        try
        {
            VMAParser parser = new VMAParser();
            DPGPAlignment ali = parser.parseFiles(vma, iCutoff, specs, strains);
            if(ali==null)
            {
                strErr = parser.getLastErrorString();
                lec = ErrorCode.ExecutionError;
                return null;
            }
            lec = ErrorCode.Ok;
            return ali;
        }
        catch(OutOfMemoryError oom)
        {
            lec = ErrorCode.VMError;
            strErr = "Java VM ran out of memory. Try increasing Java heap space.";
            return null;
        }
    }

    /**
     *  Returns the description of the last error occured.
     *
     *  @return
     */
    public String getLastError()
    {
        return strErr;
    }

    /**
     *  Returns the last error code.
     * 
     *  @return
     */
    public ErrorCode getLastErrorCode()
    {
        return lec;
    }

    /**
     *  Builds the dataset.
     *
     *  @param alignment
     *  @param fragments
     *  @param strChromosome
     *  @return
     */
    private Dataset buildDataset(DPGPAlignment alignment, HashMap<String, Vector<Fragment>> fragments, String strChromosome)
    {
        DatasetBuilder dsb = new DatasetBuilder();
        Dataset ds = dsb.build(alignment, fragments, strChromosome);
        if(ds==null)
        {
            strErr = dsb.getLastErrorString();
            lec = ErrorCode.ExecutionError;
        }
        id.wd.close();
        lec = ErrorCode.Ok;
        return ds;
    }
}
