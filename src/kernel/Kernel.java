/*
    File:
        Kernel.java
 *   
    Revision:
        2.3.0.1
 * 
    Description:
        Application kernel.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package kernel;

import algorithms.alignment.GeneFusioner;
import algorithms.alignment.GotohPairwiseFreeshiftAligner;
import dpgp.DPGPImporter;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.dna.CustomCodonTable;
import bio.gene.dna.DefaultCodonTable;
import bio.gene.dna.ICodonTable;
import gui.IWaitDialog;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Locale;
import java.util.Vector;
import plugin.loader.PluginLoader;
import plugin.PluginType;
import plugin.classes.ADatasetAnalyzer;
import plugin.classes.ADatasetExporter;
import plugin.classes.ADatasetFilter;
import plugin.classes.ADatasetImporter;
import plugin.classes.IAligner;
import plugin.classes.IGAPlugin;
import plugin.loader.AlignerLoader;


public class Kernel
{
    private final String DIR_IMPORTERS = "Plugins"+File.separator+"Importers"+File.separator;
    private final String DIR_EXPORTERS = "Plugins"+File.separator+"Exporters"+File.separator;
    private final String DIR_ANALYZERS = "Plugins"+File.separator+"Analyzers"+File.separator;
    private final String DIR_FILTERS   = "Plugins"+File.separator+"Filters"+File.separator;
    private final String DIR_ALIGNERS  = "Plugins"+File.separator+"Aligners"+File.separator;
    private final String DIR_TABLES    = "Codon tables"+File.separator;
    private final String SETTINGSFILE  = "Settings.xml";

    private Vector<ADatasetImporter>     importers      = null;
    private Vector<ADatasetExporter>     exporters      = null;
    private Vector<ADatasetFilter>       filters        = null;
    private Vector<ADatasetAnalyzer>     analyzers      = null;
    private Vector<IAligner>             aligners       = null;
    
    private Vector<ICodonTable>          codontables    = null;
    

    private QualityChecker qc                           = null;
    
    private InitData initData                           = null;
    private int iInitTable                              = -1;   // Index of the initial codon table

    private Dataset ds                                  = null;
    private int[] selInd                                = null;

    private String strLastError                         = null;

    
    public ErrorCode initialize()
    {
        // Initialize vectors.
        importers   = new Vector<ADatasetImporter>();
        exporters   = new Vector<ADatasetExporter>();
        filters     = new Vector<ADatasetFilter>();
        analyzers   = new Vector<ADatasetAnalyzer>();
        codontables = new Vector<ICodonTable>();
        aligners    = new Vector<IAligner>();

        // Load codon tables.
        if(loadCodonTables(DIR_TABLES)==ErrorCode.DirectoryDoesNotExist)
            (new File(DIR_TABLES)).mkdirs();

        // Initialization data.
        initData = new InitData(this);
        initData.sm = SettingsManager.create(SETTINGSFILE);
        // Since the kernel is initialized before the GUI, all calls to the wait dialog
        // will fail. Thus, create a dummy waiting dialog, which intercepts this calls.
        initData.wd = new IWaitDialog() {
            public void show(TYPE type) {}
            public void setText(String strMainText, String strHintText) {}
            public void close() {}
        };

        initData.ct = initializeCodonTable(initData.sm.getSetting("", SettingsManager.CODONTABLE));
        initData.locale = initializeLocale(initData.sm.getSetting("", SettingsManager.LOCALE));
        initData.codontables = listCodonTables();
        initData.aligners = listAligners();
                
        qc = new QualityChecker(initData);

        // Load plugins.
        AlignerLoader al = new AlignerLoader();        
        if(loadAligners(DIR_ALIGNERS, al)==ErrorCode.DirectoryDoesNotExist)
            (new File(DIR_ALIGNERS)).mkdirs();
        for(IAligner ali:aligners)
            ali.initialize(initData);
        
        PluginLoader pl = new PluginLoader();
        if(loadPlugins(DIR_IMPORTERS, pl)==ErrorCode.DirectoryDoesNotExist)
            (new File(DIR_IMPORTERS)).mkdirs();
        if(loadPlugins(DIR_EXPORTERS, pl)==ErrorCode.DirectoryDoesNotExist)
            (new File(DIR_EXPORTERS)).mkdirs();
        if(loadPlugins(DIR_FILTERS, pl)==ErrorCode.DirectoryDoesNotExist)
            (new File(DIR_FILTERS)).mkdirs();
        if(loadPlugins(DIR_ANALYZERS, pl)==ErrorCode.DirectoryDoesNotExist)
            (new File(DIR_ANALYZERS)).mkdirs();
        
        // Readers.
        importers.add(new builtin.readers.nfa.PluginMain());
        importers.add(new builtin.readers.jpf.PluginMain());
        importers.add(new builtin.readers.paf.PluginMain());
        importers.add(new builtin.readers.anc.PluginMain());
        for(ADatasetImporter dsi:importers)
            dsi.Initialize(initData);

        // Writers.
        exporters.add(new builtin.writers.nfa.PluginMain());
        exporters.add(new builtin.writers.nfa2.PluginMain());
        exporters.add(new builtin.writers.paf.PluginMain());
        exporters.add(new builtin.writers.gre.PluginMain());
        exporters.add(new builtin.writers.greex.PluginMain());
        exporters.add(new builtin.writers.cgs.PluginMain());
        exporters.add(new builtin.writers.embl.PluginMain());        
        for(ADatasetExporter dse:exporters)
            dse.Initialize(initData);

        // Filters.
        filters.add(new builtin.filters.shortIntrons.PluginMain());
        filters.add(new builtin.filters.chromosome.PluginMain());
        filters.add(new builtin.filters.name.PluginMain());
        filters.add(new builtin.filters.overlap.PluginMain());
        filters.add(new builtin.filters.quality.PluginMain());
        filters.add(new builtin.filters.empty.PluginMain());
        filters.add(new builtin.filters.splicing.PluginMain());
        for(ADatasetFilter dsf:filters)
            dsf.Initialize(initData);

        // Analyses.
        analyzers.add(new builtin.analyses.ffd.PluginMain());
        analyzers.add(new builtin.analyses.introns.PluginMain());
        analyzers.add(new builtin.analyses.synnonsyn.PluginMain());
        analyzers.add(new builtin.analyses.daf.PluginMain());
        analyzers.add(new builtin.analyses.composition.PluginMain());
        analyzers.add(new builtin.analyses.subst.PluginMain());
        analyzers.add(new builtin.analyses.indel.PluginMain());
        for(ADatasetAnalyzer dsa:analyzers)
            dsa.Initialize(initData);

        // Aligners.
        aligners.add(new GotohPairwiseFreeshiftAligner());
        for(IAligner aligner:aligners)
            aligner.initialize(initData);

        return ErrorCode.Ok;
    }

    /**
     *  Updates the settings after they were changed. This method should be called
     *  by the MainForm only!
     */
    public void updateSettings()
    {
        initData.locale = initializeLocale(initData.sm.getSetting("", SettingsManager.LOCALE));
        initData.ct = initializeCodonTable(initData.sm.getSetting("", SettingsManager.CODONTABLE));
    }

    /************************************************************************
    *                   DATASET OPERATIONS.                                 *
    ************************************************************************/
    /**
     *  Loads the files using the specified loader. If multiple files are specified
     *  the method combines them to one data set.
     *  If no dataset importer with the specified name is loaded, the method returns
     *  ErrorCode.ObjectNotFound.
     *
     *  @param files
     *  @param strLoader
     *  @param strParams
     *  @return
     */
    public ErrorCode importDataset(File[] files, String strLoader, String strParams)
    {
        // Search for the loader name.
        for(int i=0;i<importers.size();i++)
        {
            if(importers.get(i).GetName().equalsIgnoreCase(strLoader))
                return importDataset(files, i, strParams);
        }
        return ErrorCode.ObjectNotFound;
    }

    /**
     *  Imports DPGP dataset.
     *
     *  Note:
     *  If a dataset is loaded, it is unloaded, before the data is imported.
     *
     *  @param strParams
     *  @return
     */
    public ErrorCode importDPGPDataset(String strParams)
    {
        DPGPImporter imp = new DPGPImporter();
        if(imp.initialize(initData)!=ErrorCode.Ok)
        {
            strLastError = "DPGP Importer could not be initialized";
            return ErrorCode.ExecutionError;
        }
        Dataset tmp = imp.importDataset(strParams);
        if(tmp==null)
        {
            strLastError = imp.getLastError();
            return imp.getLastErrorCode();
        }
        initData.wd.show(IWaitDialog.TYPE.Kernel);
        ds = tmp;
        qc.validateDataset(ds);
        selectAll();
        initData.wd.close();
        return ErrorCode.Ok;
    }

    public ErrorCode fuseGeneEntries(int[] indices, String strParams)
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        // Create new data set.
        Dataset tmp = new Dataset();
        if(indices==null)
            indices = selInd;
        for(int i:indices)
        {
            if(i>-1 && i<ds.getGenesCount())
                tmp.addGene(ds.getGeneEntry(i));
        }
        if(tmp.getGenesCount()==0)
        {
            strLastError = "No gene entries selected to fuse";
            return ErrorCode.SelectionIsEmpty;
        }
        // Fuse the genes.
        GeneFusioner gf = new GeneFusioner(initData.wd);
        tmp = gf.fuseGenes(tmp, listAligners());
        if(tmp==null)
            return ErrorCode.CancelledByUser;
        else if(tmp.getGenesCount()>0)
        {
            int nOld = ds.getGenesCount();
            ds.merge(tmp);
            qc.validateDataset(ds);
            if(selInd==null)
            selectAll();
            else
            {
                int[] ind = Arrays.copyOf(selInd, selInd.length+ds.getGenesCount()-nOld);
                for(int i=selInd.length;i<ind.length;i++)
                    ind[i] = nOld++;
                selInd = ind;
            }
        }
        return ErrorCode.Ok;
    }

    /**
     *  Adds the new gene entry to the loaded dataset. If no dataset is loaded,
     *  an empty dataset is created and the entry is added to it. If the dataset
     *  already has the gene with the same common name and bForce is false,
     *  the method returns ErrorCode.ObjectExists.
     *
     *  @param ge
     *  @return
     */
    public ErrorCode addGeneEntry(GeneEntry ge, boolean bForce)
    {
        if(ds==null)
            ds = new Dataset();
        if(ds.hasGene(ge.getCommonName()) && !bForce)
            return ErrorCode.ObjectExists;
        else
        {
            ds.addGene(ge);
            qc.validateDataset(ds);
            return ErrorCode.Ok;
        }
    }

    /**
     *  Unloads the dataset.
     */
    public void unloadDataset()
    {
        ds = null;
        selInd = null;
    }

    /**
     *  Exports the specified entries of the data set into the file. If indices
     *  is null, the currently selected genes are saved. If the value in the indices
     *  array is invalid, it is ignored.
     *  If no data set exporter with the specified name is loaded, the method returns
     *  ErrorCode.ObjectNotFound.
     *
     *  Remarks:
     *      Currently selected genes are either all genes or the genes selected
     *      by the previously used filter.
     *
     *  @param indices
     *  @param iIndex
     *  @param file
     *  @param strParams
     *  @return
     */
    public ErrorCode exportDataset(int[] indices, String strName, File file, String strParams)
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        for(int i=0;i<exporters.size();i++)
        {
            if(exporters.get(i).GetName().equalsIgnoreCase(strName))
                return exportDataset(indices, i, file, strParams);
        }
        return ErrorCode.ObjectNotFound;
    }    

    /**
     *  Applies the specified filter to the data set.
     *  If no dataset filter with the specified name is loaded, the method
     *  returns ErrroCode.ObjectNotFound.
     *
     *  @param strName
     *  @param strParams
     *  @return
     */
    public ErrorCode applyFilter(String strName, String strParams)
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        for(int i=0;i<filters.size();i++)
        {
            if(filters.get(i).GetName().equalsIgnoreCase(strName))
                return applyFilter(i, strParams);
        }
        return ErrorCode.ObjectNotFound;
    }    

    /**
     *  Analyzes the specified genes with the dataset analyzer specified by its
     *  name. If indices is null, currently selected genes are analyzed.
     *  
     *  Remarks:
     *      Currently selected genes are either all genes or the genes selected
     *      by the previously used filter.
     * 
     *  @param indices
     *  @param strName
     *  @param strParams
     *  @return
     */
    public ErrorCode performAnalysis(int[] indices, String strName, String strParams)
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        for(int i=0;i<analyzers.size();i++)
        {
            if(analyzers.get(i).GetName().equalsIgnoreCase(strName))
                return performAnalysis(indices, i, strParams);
        }
        return ErrorCode.ObjectNotFound;
    }    

    /**
     *  Returns currently selected indices.
     * 
     *  Note:
     *  If no dataset is loaded, the method returns null.
     *
     *  @return
     */
    public int[] getSelectedIndices()
    {
        return selInd;
    }

    /**
     *  Selects all gene entries in the loaded dataset.
     */
    public void selectAll()
    {
        if(ds==null)
            return;
        selInd = new int[ds.getGenesCount()];
        for(int i=0;i<ds.getGenesCount();i++)
            selInd[i] = i;
    }
    
    /**
     *  Inverts the current selection.
     */
    public void invertSelection()
    {
        if(ds==null)
            return;
        int[] tmp = new int[ds.getGenesCount()];
        for(int i=0;i<tmp.length;i++)
            tmp[i] = i;
        for(int i:selInd)
            tmp[i]=-1;
        int[] result = new int[tmp.length-selInd.length];
        int iIndex = 0;
        for(int i:tmp)
        {
            if(i>-1)
            {
                result[iIndex] = i;
                iIndex++;
            }
        }
        selInd = result;
    }

    /**
     *  Replaces the current selection.
     *
     *  @param indices
     */
    public ErrorCode setSelectedIndices(int[] indices)
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        selInd = indices;
        return ErrorCode.Ok;
    }    
    
    /** 
     *  Returns the loaded dataset.
     * 
     *  @return
     */
    public Dataset getDataset()
    {
        return ds;
    }    
    
    /**
     *  Returns true if the data set is loaded.
     * 
     *  @return
     */
    public boolean isDatasetLoaded()
    {
        return ds!=null;
    }

    /**
     *  Sorts the genes of the dataset by name.
     *
     *  NOTE: after sort the selection is lost!
     */
    public void sort()
    {
        if(ds==null)
            return;
        ds.sort();
        selectAll();
    }

    /**
     *  Sorts the genes of the dataset by the quality level.
     *
     *  NOTE: after sort the selection is lost!
     */
    public void sortByQualityLevel()
    {
        if(ds==null)
            return;
        ds.sortByQuality();
        selectAll();
    }

    /**
     *  Runs quality check on the loaded dataset, if any.
     *
     *  @return
     */
    public ErrorCode runQualityCheck()
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        qc.validateDataset(ds);
        return ErrorCode.Ok;
    }

    /************************************************************************
    *                   CODON TABLES OPERATIONS.                            *
    ************************************************************************/
    
    /**
     *  Adds the codon table to the codon tables list.
     * 
     *  Returns:
     *      Ok                  if the operation succeeds
     *      InvalidParameter    if ct is null
     *      ObjectExists        if the codon table with the same name exists
     * 
     *  @param ct
     *  @return
     */
    public ErrorCode addCodonTable(CustomCodonTable ct)
    {
        if(ct==null)
            return ErrorCode.InvalidParameter;
        // Check whether a codon table for the same organism is already loaded.
        for(ICodonTable table:codontables)
        {
            if(table.getName().equalsIgnoreCase(ct.getName()))
                return ErrorCode.ObjectExists;
        }
        codontables.add(ct);
        initData.codontables = listCodonTables();
        ct.saveToFile(DIR_TABLES+ct.getName()+CustomCodonTable.CODONTABLE_EXTENTION);
        return ErrorCode.Ok;
    }
    
    /**
     *  Returns the names of the loaded codons tables.
     *
     *  @return
     */
    public String[] listCodonTables()
    {
        String[] names = new String[codontables.size()];
        for(int i=0;i<codontables.size();i++)
            names[i] = codontables.get(i).getName();
        return names;
    }

    /**
     *  Searches the loaded codon tables for the codon table with the specified
     *  name, and, if succeeds, sets the current codon table to the found one.
     *  The initial codon table can be restored by providing "Initial" in strName.
     *  Providing "Default" in strName will set the coding table to the
     *  universal codon table.
     *
     *  @param strName
     *  @return
     */
    public ErrorCode setCodonTable(String strName)
    {
        if(strName.equalsIgnoreCase("Initial"))
        {
            initData.ct = codontables.get(iInitTable);
            return ErrorCode.Ok;
        }
        else if(strName.equalsIgnoreCase("Default"))
        {
            initData.ct = codontables.get(0);
            return ErrorCode.Ok;
        }
        else
        {
            for(ICodonTable ct:codontables)
            {
                if(ct.getName().equalsIgnoreCase(strName))
                {
                    initData.ct = ct;
                    return ErrorCode.Ok;
                }
            }
        }
        return ErrorCode.ObjectNotFound;
    }

    /************************************************************************
    *                      PLUGINS OPERATIONS.                              *
    ************************************************************************/
    /**
     *  Returns the array of PluginDetails objects, containing the details
     *  about the loaded plugins of the specified type.
     *
     *  @param type
     *  @return
     */
    public PluginDetails[] getPluginsDetails(PluginType type)
    {
        // Declare the array.
        PluginDetails[] details = null;
        // Vector to iterate through.
        Vector v = null;
        // Distinguish the plugin type.
        switch(type)
        {
            case ANALYZER: v = analyzers; break;
            case EXPORTER: v = exporters; break;
            case IMPORTER: v = importers; break;
            case FILTER  : v = filters; break;
        }
        // Iterate through the vector.
        details = new PluginDetails[v.size()];
        for(int i=0;i<v.size();i++)
        {
            details[i] = new PluginDetails();
            details[i].strMenuItemName = ((IGAPlugin)(v.get(i))).GetMenuItemName();
            details[i].strName = ((IGAPlugin)(v.get(i))).GetName();
            details[i].strDescription = ((IGAPlugin)(v.get(i))).GetDescription();
            details[i].strParamStr = ((IGAPlugin)(v.get(i))).GetParamString();
            details[i].bSmd = ((IGAPlugin)(v.get(i))).SupportsMissingData();
            details[i].bSad = ((IGAPlugin)(v.get(i))).SupportsAmbiguousData();
            if(type==PluginType.EXPORTER)
            {
                details[i].ed = new ExtensionDetails(exporters.get(i).GetFileExtension(),
                                                     exporters.get(i).GetFileDescription(), i);
            }
            else if(type==PluginType.IMPORTER)
            {
                details[i].ed = new ExtensionDetails(importers.get(i).GetFileExtension(),
                                                     importers.get(i).GetFileDescription(), i);
            }
            else
                details[i].ed = null;
        }
        return details;
    }

    /**
     *  Lists all loaded aligners.
     * 
     *  @return
     */
    public IAligner[] listAligners()
    {
        return aligners.toArray(new IAligner[1]);
    }

    /**
     *  Returns the initialization data.
     *
     *  @return
     */
    public InitData getInitializationData()
    {
        return initData;
    }

    /**
     *  Sets the wait dialog.
     *
     *  @param wd
     *  @return
     */
    public ErrorCode setWaitDialog(IWaitDialog wd)
    {
        if(wd==null)
            return ErrorCode.InvalidParameter;
        initData.wd = wd;
        return ErrorCode.Ok;
    }

    /**
     *  Returns the string description of the last error occured.
     *
     *  @return
     */
    public String getLastErrorString()
    {
        return strLastError;
    }

    /************************************************************************
    *                   PRIVATE METHODS.                                    *
    ************************************************************************/
    private ErrorCode loadPlugins(String strDir, PluginLoader pl)
    {
        File dir = new File(strDir);
        // Check whether the file is a valid existing directory.
        if(!dir.isDirectory() || !dir.exists())
            return ErrorCode.DirectoryDoesNotExist;
        String[] children = dir.list(new FilenameFilter()
         {
            public boolean accept(File dir, String name)
            {
                return name.endsWith(IGAPlugin.PLUGIN_EXTENTION);
            }
         });
        for(String s:children)
        {
            IGAPlugin plugin = pl.loadPlugin(dir.getAbsolutePath()+File.separator+s);
            plugin.Initialize(initData);
            // Determine the type of the plugin.
            switch(plugin.GetType())
            {
                case IMPORTER: importers.add((ADatasetImporter)plugin); break;
                case EXPORTER: exporters.add((ADatasetExporter)plugin); break;
                case ANALYZER: analyzers.add((ADatasetAnalyzer)plugin); break;
                case FILTER  : filters.add((ADatasetFilter)plugin); break;
            }
        }
        return ErrorCode.Ok;
    }

    /**
     *  Loads the codon tables.
     *
     *  @param strDir
     *  @return
     */
    public ErrorCode loadCodonTables(String strDir)
    {
        // First, add the built-in codon table.
        codontables.add(new DefaultCodonTable());
        File dir = new File(strDir);
        // Check whether the file is a valid existing directory.
        if(!dir.isDirectory() || !dir.exists())
            return ErrorCode.DirectoryDoesNotExist;
        String[] children = dir.list(new FilenameFilter()
                 {
                    public boolean accept(File dir, String name)
                    {
                        return name.endsWith(CustomCodonTable.CODONTABLE_EXTENTION);
                    }
                 });
        for(String s:children)
        {
            CustomCodonTable cct = CustomCodonTable.loadFromFile(strDir+File.separator+s);
            if(cct!=null)
                codontables.add(cct);
        }
        return ErrorCode.Ok;
    }

    /**
     *  Loads the aligners.
     *
     *  @param strDir
     *  @return
     */
    private ErrorCode loadAligners(String strDir, AlignerLoader al)
    {
        File dir = new File(strDir);
        // Check whether the file is a valid existing directory.
        if(!dir.isDirectory() || !dir.exists())
            return ErrorCode.DirectoryDoesNotExist;
        String[] children = dir.list(new FilenameFilter()
         {
            public boolean accept(File dir, String name)
            {
                return name.endsWith(IAligner.PLUGIN_EXTENTION);
            }
         });
        for(String s:children)
        {
            IAligner aligner = al.loadPlugin(dir.getAbsolutePath()+File.separator+s);
            aligners.add(aligner);
        }
        return ErrorCode.Ok;
    }

    /**
     *  Initializes the initial codon table. If the codon table with the specified
     *  name cannot be found, default codon table is used instead.
     *
     *  @param strName
     *  @return
     */
    private ICodonTable initializeCodonTable(String strName)
    {
        if(strName==null || strName.isEmpty())
        {
            iInitTable = 0;
            return codontables.get(0);
        }
        // Try to locate the specified codon table in the list of
        // loaded codon tables.
        for(int i=0;i<codontables.size();i++)
        {
            ICodonTable ct = codontables.get(i);
            if(ct.getName().equalsIgnoreCase(strName))
            {
                iInitTable = i;
                return ct;
            }
        }
        // If the codon table could not be found, initialize the value with the
        // default codon table.
        iInitTable = 0;
        return codontables.get(0);
    }

    private Locale initializeLocale(String strCountry)
    {
        if(strCountry==null || strCountry.isEmpty())
            return Locale.GERMANY;
        if(strCountry.equals("Germany"))
            return Locale.GERMANY;
        else
            return Locale.US;
    }

    /**
     *  Loads the files using the specified loader. If multiple files are specified
     *  the method combines them to one data set. The loader is specified by its
     *  kernel index.
     *  If the index value is invalid, the method returns ErrorCode.ObjectNotFound.
     *
     *  @param files
     *  @param strLoader
     *  @param strParams
     *  @return
     */
    private ErrorCode importDataset(File[] files, int iIndex, String strParams)
    {
        if(iIndex>=importers.size() || iIndex<0)
            return ErrorCode.ObjectNotFound;
        // Files are assumed to be existent and of the type the specified
        // plugin can load.
        if(ds==null)
            ds = new Dataset();
        Dataset tmp = null;
        try
        {
            tmp = importers.get(iIndex).ImportDataset(files, strParams);
            if(tmp==null)
            {
                strLastError = importers.get(iIndex).GetLastError();
                return ErrorCode.ExecutionError;
            }
        }
        catch(Exception e)
        {
            strLastError = e.getMessage();
            initData.wd.close();
            return ErrorCode.ExecutionError;
        }
        int nOld = ds.getGenesCount();
        initData.wd.show(IWaitDialog.TYPE.Kernel);
        ds.merge(tmp);        
        qc.validateDataset(ds);
        initData.wd.close();
        if(selInd==null)
            selectAll();
        else
        {
            int[] ind = Arrays.copyOf(selInd, selInd.length+ds.getGenesCount()-nOld);
            for(int i=selInd.length;i<ind.length;i++)
                ind[i] = nOld++;
            selInd = ind;
        }
        return ErrorCode.Ok;
    }

    /**
     *  Saves the specified entries of the data set into the file. If indices
     *  is null, the currently selected genes are saved. If the value in the indices
     *  array is invalid, it is ignored.
     *  If the index value is invalid, the method returns ErrorCode.ObjectNotFound.
     *
     *  Remarks:
     *      Currently selected genes are either all genes or the genes selected
     *      by the previously used filter.
     *
     *  @param indices
     *  @param iIndex
     *  @param file
     *  @param strParams
     *  @return
     */
    private ErrorCode exportDataset(int[] indices, int iIndex, File file, String strParams)
    {
        if(iIndex>=exporters.size() || iIndex<0)
            return ErrorCode.ObjectNotFound;
        // Create new data set.
        Dataset tmp = new Dataset();
        if(indices==null)
            indices = selInd;
        for(int i:indices)
        {
            if(i>-1 && i<ds.getGenesCount())
                tmp.addGene(ds.getGeneEntry(i));
        }
        if(tmp.getGenesCount()==0)
        {
            strLastError = "No valid data to export";
            return ErrorCode.SelectionIsEmpty;
        }
        // Save the dataset.
        try
        {
            ErrorCode ec = exporters.get(iIndex).ExportDataset(tmp, file, strParams);
            if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
            {
                strLastError = exporters.get(iIndex).GetLastError();
                return ErrorCode.ExecutionError;
            }
        }
        catch(Exception e)
        {
            strLastError = e.getMessage();
            initData.wd.close();
            return ErrorCode.ExecutionError;
        }
        return ErrorCode.Ok;
    }

    /**
     *  Applies the filter specified by its index.
     *  If the index is invalid, the method returns ErrorCode.ObjectNotFound.
     *
     *  Remarks:
     *      If the dataset was filtered before, only the selected indices are
     *      used by the current filter. Thus, the filters can be combined.
     *
     *  @param iIndex
     *  @param strParams
     *  @return
     */
    private ErrorCode applyFilter(int iIndex, String strParams)
    {
        if(ds==null)
            return ErrorCode.NoDatasetLoaded;
        if(iIndex>=filters.size() || iIndex<0)
            return ErrorCode.ObjectNotFound;
        try
        {
            int[] tmp = filters.get(iIndex).GetFilteredIndices(ds, strParams);
            if(tmp==null)
            {
                strLastError = filters.get(iIndex).GetLastError();
                return ErrorCode.ExecutionError;
            }
            if(tmp.length==1 && tmp[0]==-1)
                return ErrorCode.CancelledByUser;
            // Combine the filters.
            Vector<Integer> ind = new Vector<Integer>();
            for(int i:selInd)
            {
                for(int j:tmp)
                {
                    if(i==j)
                        ind.add(i);
                }
            }
            selInd = new int[ind.size()];
            for(int i=0;i<selInd.length;i++)
                selInd[i] = ind.get(i);
        }
        catch(Exception e)
        {
            strLastError = e.getMessage();
            initData.wd.close();
            return ErrorCode.ExecutionError;
        }
        return ErrorCode.Ok;
    }

    /**
     *  Performs the specified analysis on the specified genes.
     *  If indices is null, currently selected genes are used.
     *
     *  Remarks:
     *      Currently selected genes are either all genes or the genes selected
     *      by the previously used filter.
     *
     *  @param indices
     *  @param iIndex
     *  @param strParams
     *  @return
     */
    private ErrorCode performAnalysis(int[] indices, int iIndex, String strParams)
    {
        if(iIndex>=analyzers.size() || iIndex<0)
            return ErrorCode.ObjectNotFound;
        // Create new data set.
        Dataset tmp = new Dataset();
        if(indices==null)
            indices = selInd;
        for(int i:indices)
        {
            if(i>-1 && i<ds.getGenesCount())
                tmp.addGene(ds.getGeneEntry(i));
        }
        if(tmp.getGenesCount()==0)
        {
            strLastError = "No valid data to analyze";
            return ErrorCode.SelectionIsEmpty;
        }
        try
        {
            ErrorCode ec = analyzers.get(iIndex).AnalyzeDataset(tmp, strParams);
            if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
            {
                strLastError = analyzers.get(iIndex).GetLastError();
                return ErrorCode.ExecutionError;
            }
        }
        catch(Exception e)
        {
            strLastError = e.getMessage();
            initData.wd.close();
            return ErrorCode.ExecutionError;
        }
        return ErrorCode.Ok;
    }
}
