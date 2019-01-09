/*
    File:
        MainForm.java
 *   
    Revision:
        2.6.0.1
 * 
    Description:
        Application main form.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui;

import algorithms.Coverage;
import gui.gene.AddGeneDialog;
import gui.script.ScriptForm;
import javax.swing.event.ChangeEvent;
import kernel.PluginDetails;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.dna.CustomCodonTable;
import gui.editor.EditorForm;
import gui.views.GeneAnnotationView;
import gui.views.GeneEntryDetailsView;
import java.awt.Cursor;
import plugin.classes.AGeneAnalyzerView;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.FilenameFilter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;
import javax.swing.ButtonGroup;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTable;
import javax.swing.JToggleButton;
import javax.swing.event.ChangeListener;
import javax.swing.table.DefaultTableModel;
import kernel.fix.DatasetFixer;
import kernel.ErrorCode;
import kernel.ExtensionDetails;
import kernel.InitData;
import kernel.Kernel;
import kernel.Main;
import kernel.QualityChecker;
import kernel.script.ScriptInterpreter;
import kernel.SettingsManager;
import plugin.PluginType;
import plugin.classes.IAligner;
import plugin.loader.ViewLoader;


public class MainForm extends JFrame
{
    // Monitors the click event on a menu item and calls the routine depending on
    // the type of plugin selected.
    private class PluginMenuItemListener implements ActionListener
    {
        // Index of the plugin to use.
        private int             index;
        private PluginType      type;
        
        public PluginMenuItemListener(int index, PluginType type)
        {
            this.index = index;
            this.type  = type;
        }
        
        public void actionPerformed(ActionEvent e)
        {
            switch(type)
            {
                case ANALYZER: performAnalysis(index); break;
                case FILTER  : applyFilter(index); break;
                default : return;   // Any other plugin type.
            }
        }
    };
    
    private class CTMenuItemListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
            kernel.setCodonTable(e.getActionCommand());
        }
    };


    private final String DIR_VIEWS              = "Plugins"+File.separator+"Views"+File.separator;

    private static final String WIN_X           = "WindowX";
    private static final String WIN_Y           = "WindowY";
    private static final String WIN_STATE       = "WindowState";
    private static final String WIN_WIDTH       = "WindowWidth";
    private static final String WIN_HEIGHT      = "WindowHeight";

    
    private Kernel kernel                       = null;
    private ScriptForm sf                       = null;
    private ScriptInterpreter interpreter       = null;
    private CodonTableCreator ctc               = null;
    private SettingsManager sm                  = null;
    private InitData  initData                  = null;
    private DatasetFixer dsf                    = null;

    // Views.
    private AGeneAnalyzerView[] views           = null;

    // Plugins.
    private PluginDetails[] importers           = null;
    private PluginDetails[] exporters           = null;
    private PluginDetails[] analyzers           = null;
    private PluginDetails[] filters             = null;

    /**
     *  Creates the main form and initializes the controls.
     */
    public MainForm(Kernel kernel) 
    {
        // Set application title.
        setTitle(Main.APPTITLE);
        Thread.currentThread().setName("GeneAnalyzer Main Form");
        // Create kernel.
        this.kernel = kernel; 
        this.interpreter = new ScriptInterpreter(kernel, this);
        initData = kernel.getInitializationData();
        sm = initData.sm;
        this.setGlassPane(new WaitDialog());
        kernel.setWaitDialog((IWaitDialog)this.getGlassPane());
        // Initialize controls.
        initComponents();
        setLocationRelativeTo(null);
        setMinimumSize(new Dimension(1100, 700));
             
        // Initialize list view.
        lvGenes.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        lvGenes.getColumnModel().getColumn(0).setPreferredWidth(40);   // Use
        lvGenes.getColumnModel().getColumn(1).setPreferredWidth(40);   // #
        lvGenes.getColumnModel().getColumn(2).setPreferredWidth(95);   // Common name
        lvGenes.getColumnModel().getColumn(3).setPreferredWidth(95);   // Alias
        lvGenes.getColumnModel().getColumn(4).setPreferredWidth(50);   // Strains count
        lvGenes.getColumnModel().getColumn(5).setPreferredWidth(132);  // Quality check

        for(PluginType type:PluginType.values())
        {
            JMenu parent = null;
            switch(type)
            {
                // Skip importers and exporters.
                case IMPORTER:
                    importers = kernel.getPluginsDetails(type);
                    continue;
                case EXPORTER: 
                    exporters = kernel.getPluginsDetails(type);
                    continue;
                // Update analyzers and filters menus.
                case ANALYZER:
                    analyzers = kernel.getPluginsDetails(type);
                    parent = jmAnalyses;
                    break;
                case FILTER  :
                    filters = kernel.getPluginsDetails(type);
                    parent = jmFilters;
                    break;
            }
            // Add items.
            PluginDetails[] plugins = kernel.getPluginsDetails(type);
            for(int i=0;i<plugins.length;i++)
            {
                JMenuItem mi = new JMenuItem(plugins[i].strMenuItemName);
                mi.setToolTipText(plugins[i].strDescription);
                mi.addActionListener(new PluginMenuItemListener(i, type));
                parent.add(mi);
            }
        }
        // Initialize views.
        views = loadViews(DIR_VIEWS);
        // Add views.
        ButtonGroup bg = new ButtonGroup();
        for(int i=0;i<views.length;i++)
        {
            jpToolView.add(views[i]);
            views[i].setVisible(false);
            final JToggleButton jtb = new JToggleButton();
            jtb.setText(views[i].getButtonName());
            jtb.setToolTipText(views[i].getButtonHintText());
            final int index = i;
            jtb.addChangeListener(new ChangeListener() {
                public void stateChanged(ChangeEvent e)
                {
                    views[index].setVisible(jtb.isSelected());
                }
                }); 
            bg.add(jtb);
            jViewsToolbar.add(jtb);
            if(i==0)
                jtb.setSelected(true);
        }
        views[0].setVisible(true);
        jpToolView.validate();        
        // Set window position.       
        String x = (String)sm.getSetting("", WIN_X);
        String y = (String)sm.getSetting("", WIN_Y);
        if(x!=null && y!=null)
        {
            Point p = new Point(Integer.parseInt(x), Integer.parseInt(y));
            this.setLocation(p);
        }
        // Set window size.
        String w = (String)sm.getSetting("", WIN_WIDTH);
        String h = (String)sm.getSetting("", WIN_HEIGHT);
        if(w!=null && h!=null)
        {
            Dimension d = new Dimension(Integer.parseInt(w), Integer.parseInt(h));
            this.setSize(d);
        }
        // Set window state.
        String state = (String)sm.getSetting("", WIN_STATE);
        if(state!=null)
            this.setExtendedState(Integer.parseInt(state));
    }

    /**
     *  Updates the gene list.
     */
    public void updateGeneList()
    {
        DefaultTableModel model = (DefaultTableModel)lvGenes.getModel();
        model.setRowCount(0);
        // Add new genes.
        Dataset ds = kernel.getDataset();
        if(ds==null)
            return;
        for(int i=0;i<ds.getGenesCount();i++)
        {
            String strQuality = null;
            Integer quality = (Integer)ds.getGeneEntry(i).getProperty(QualityChecker.QUALITY_LEVEL);
            switch(quality)
            {
                case QualityChecker.QUALITY_LEVEL_0:
                    strQuality = "OK";
                    break;
                case QualityChecker.QUALITY_LEVEL_1:
                    strQuality = "Negligible errors";
                    break;
                case QualityChecker.QUALITY_LEVEL_2:
                    strQuality = "Invalid ORF";
                    break;
                case QualityChecker.QUALITY_LEVEL_3:
                    strQuality = "Premature terminal codon";
                    break;
                case QualityChecker.QUALITY_LEVEL_4:
                    strQuality = "Severe deficits";
                    break;
                case QualityChecker.QUALITY_LEVEL_5:
                    strQuality = "Fatal errors";
                    break;
            }
            model.addRow(new Object[]{new Boolean(false),
                         i+1,
                         ds.getGeneEntry(i).getCommonName(),
                         ds.getGeneEntry(i).getAlias(),
                         ds.getGeneEntry(i).getStrainsCount(),
                         strQuality});
        }
        int[] tmp = kernel.getSelectedIndices();
        if(tmp==null)
            return;
        for(int i:tmp)
            model.setValueAt(true, i, 0);
    }    

    /**
     *  Applies the filter with specified index to the loaded data set.
     *  If indices is not null, index is ignored and the indices are used.
     * 
     *  @param index
     *  @param indices
     */
    private void applyFilter(final int index)
    {
        Runnable r = new Runnable()
        {
            public void run()
            {
                // Filter the data.
                kernel.setSelectedIndices(getSelectedIndices());
                int[] indices = null;
                ErrorCode ec = kernel.applyFilter(filters[index].strName, null);
                if(ec==ErrorCode.Ok)
                   indices = kernel.getSelectedIndices();
                else if(ec!=ErrorCode.CancelledByUser)
                {
                    String strErrDesc = kernel.getLastErrorString();
                    String strErrMsg = null;
                    if(strErrDesc==null || strErrDesc.isEmpty())
                        strErrMsg = String.format("<html>An unspecified error occured while filtering the data with <b>%s</b></html>",
                                                  filters[index].strName);
                    else
                        strErrMsg = String.format("<html>An error occured while filtering the data with <b>%s:</b><br>\t%s</html>",
                                                  filters[index].strName, strErrDesc);
                    JOptionPane.showMessageDialog(null, strErrMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
                }
                if(indices==null)
                    return;
                DefaultTableModel model = (DefaultTableModel)lvGenes.getModel();
                // Clear all genes.
                for(int i=0;i<model.getRowCount();i++)
                    model.setValueAt(new Boolean(false), i, 0);
                // Set indices.
                for(int i=0;i<indices.length;i++)
                    model.setValueAt(new Boolean(true), indices[i], 0);
            }
        };
        new Thread(r).start();
    }
       
    private void performAnalysis(final int index)
    {
        Runnable r = new Runnable()
        {
            public void run()
            {
                ErrorCode ec = kernel.performAnalysis(getSelectedIndices(), analyzers[index].strName, null);
                if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
                {
                    String strErrDesc = kernel.getLastErrorString();
                    String strErrMsg = null;
                    if(strErrDesc==null || strErrDesc.isEmpty())
                        strErrMsg = String.format("<html>An unspecified error occured while analyzing the data with <b>%s</b></html>",
                                                  analyzers[index].strName);
                    else
                        strErrMsg = String.format("<html>An error occured while analyzing the data with <b>%s:</b><br>\t%s</html>",
                                                  analyzers[index].strName, strErrDesc);
                    JOptionPane.showMessageDialog(null, strErrMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
                }
            }
        };
        new Thread(r).start();
    }

    /**
     *  Returns the indices of selected genes.
     * 
     *  @return
     */
    public int[] getSelectedIndices()
    {
        // Indices.
        Vector<Integer> tmp = new Vector<Integer>();
        DefaultTableModel model = (DefaultTableModel)lvGenes.getModel();
        // Add the indices of genes which are selected.
        for(int i=0;i<model.getRowCount();i++)
            if((Boolean)model.getValueAt(i, 0))
                tmp.add(i);
        if(tmp.size()==0)
            return new int[]{};
        int[] indices = new int[tmp.size()];
        for(int i=0;i<tmp.size();i++)
            indices[i] = tmp.get(i);
        return indices;
    }
            
    private void updateViews()
    {
        int i = lvGenes.getSelectedRow();
        if(i>-1)
        {
            GeneEntry ge = kernel.getDataset().getGeneEntry(i);
            for(int n=0;n<views.length;n++)
                views[n].displayGeneEntryInfo(ge);
          
        }
    }

    private AGeneAnalyzerView[] loadViews(String strDir)
    {
        Vector<AGeneAnalyzerView> tmp = new Vector<AGeneAnalyzerView>();
        // Add built-in views.
        tmp.add(new GeneEntryDetailsView());
        tmp.add(new GeneAnnotationView());
        // Load plugins.
        File dir = new File(strDir);
        // Check whether the file is a valid existing directory.
        if(!dir.isDirectory() || !dir.exists())
            (new File(strDir)).mkdirs();
        ViewLoader vl = new ViewLoader();
        String[] children = dir.list(new FilenameFilter()
         {
            public boolean accept(File dir, String name)
            {
                return name.endsWith(AGeneAnalyzerView.PLUGIN_EXTENTION);
            }
         });
        for(String s:children)
        {
            AGeneAnalyzerView view = vl.loadPlugin(dir.getAbsolutePath()+File.separator+s);
            tmp.add(view);
        }
        return tmp.toArray(new AGeneAnalyzerView[tmp.size()]);
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jpmGenesMenu = new javax.swing.JPopupMenu();
        jmiSelect = new javax.swing.JMenuItem();
        jmiUnselect = new javax.swing.JMenuItem();
        jmiInvert = new javax.swing.JMenuItem();
        jSeparator6 = new javax.swing.JSeparator();
        jPanel1 = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        lvGenes = new javax.swing.JTable();
        jSeparator3 = new javax.swing.JSeparator();
        jpToolView = new javax.swing.JPanel();
        jViewsToolbar = new javax.swing.JToolBar();
        jmMain = new javax.swing.JMenuBar();
        jmDataset = new javax.swing.JMenu();
        jmImport = new javax.swing.JMenu();
        miImport = new javax.swing.JMenuItem();
        miCreateDpgp = new javax.swing.JMenuItem();
        miExport = new javax.swing.JMenuItem();
        miUnload = new javax.swing.JMenuItem();
        jSeparator5 = new javax.swing.JSeparator();
        miAddGene = new javax.swing.JMenuItem();
        miFuse = new javax.swing.JMenuItem();
        miAlign = new javax.swing.JMenuItem();
        jmFix = new javax.swing.JMenu();
        miRemoveInsertions = new javax.swing.JMenuItem();
        miRemovePTC = new javax.swing.JMenuItem();
        jSeparator4 = new javax.swing.JSeparator();
        jmCT = new javax.swing.JMenu();
        jSeparator7 = new javax.swing.JSeparator();
        miExit = new javax.swing.JMenuItem();
        jmAnalyses = new javax.swing.JMenu();
        jmFilters = new javax.swing.JMenu();
        jmTools = new javax.swing.JMenu();
        jmSort = new javax.swing.JMenu();
        miSortByName = new javax.swing.JMenuItem();
        miSortByQuality = new javax.swing.JMenuItem();
        miQualityCheck = new javax.swing.JMenuItem();
        miCreateCT = new javax.swing.JMenuItem();
        miCoverage = new javax.swing.JMenuItem();
        miShowNames = new javax.swing.JMenuItem();
        jSeparator8 = new javax.swing.JSeparator();
        miScript = new javax.swing.JMenuItem();
        jSeparator1 = new javax.swing.JSeparator();
        miSettings = new javax.swing.JMenuItem();

        jmiSelect.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.ALT_MASK | java.awt.event.InputEvent.SHIFT_MASK));
        jmiSelect.setText("Select all");
        jmiSelect.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jmiSelectActionPerformed1(evt);
            }
        });
        jpmGenesMenu.add(jmiSelect);

        jmiUnselect.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_U, java.awt.event.InputEvent.ALT_MASK | java.awt.event.InputEvent.SHIFT_MASK));
        jmiUnselect.setText("Unselect all");
        jmiUnselect.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jmiUnselectActionPerformed(evt);
            }
        });
        jpmGenesMenu.add(jmiUnselect);

        jmiInvert.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_I, java.awt.event.InputEvent.ALT_MASK | java.awt.event.InputEvent.SHIFT_MASK));
        jmiInvert.setText("Invert selection");
        jmiInvert.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jmiInvertActionPerformed(evt);
            }
        });
        jpmGenesMenu.add(jmiInvert);

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder("Data set"));

        lvGenes.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {
                "Use", "#", "Common name", "Alias", "Strains", "Quality check"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Boolean.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                true, false, false, false, false, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        lvGenes.setComponentPopupMenu(jpmGenesMenu);
        lvGenes.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        lvGenes.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                lvGenesMouseClicked(evt);
            }
        });
        lvGenes.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                lvGenesKeyReleased(evt);
            }
        });
        jScrollPane1.setViewportView(lvGenes);

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 471, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(18, Short.MAX_VALUE))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 724, Short.MAX_VALUE)
                .addContainerGap())
        );

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        jSeparator3.setOrientation(javax.swing.SwingConstants.VERTICAL);

        jpToolView.setLayout(new javax.swing.BoxLayout(jpToolView, javax.swing.BoxLayout.LINE_AXIS));

        jViewsToolbar.setFloatable(false);
        jViewsToolbar.setRollover(true);

        jmDataset.setText("Dataset");
        jmDataset.addMenuListener(new javax.swing.event.MenuListener() {
            public void menuCanceled(javax.swing.event.MenuEvent evt) {
            }
            public void menuDeselected(javax.swing.event.MenuEvent evt) {
            }
            public void menuSelected(javax.swing.event.MenuEvent evt) {
                jmDatasetMenuSelected(evt);
            }
        });

        jmImport.setText("Import");

        miImport.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_I, java.awt.event.InputEvent.ALT_MASK));
        miImport.setText("From file");
        miImport.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miImportActionPerformed(evt);
            }
        });
        jmImport.add(miImport);

        miCreateDpgp.setText("From DPGP files");
        miCreateDpgp.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCreateDpgpActionPerformed(evt);
            }
        });
        jmImport.add(miCreateDpgp);

        jmDataset.add(jmImport);

        miExport.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_E, java.awt.event.InputEvent.ALT_MASK));
        miExport.setText("Export as...");
        miExport.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miExportActionPerformed(evt);
            }
        });
        jmDataset.add(miExport);

        miUnload.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_U, java.awt.event.InputEvent.ALT_MASK));
        miUnload.setText("Unload");
        miUnload.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miUnloadActionPerformed(evt);
            }
        });
        jmDataset.add(miUnload);
        jmDataset.add(jSeparator5);

        miAddGene.setText("Add gene");
        miAddGene.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miAddGeneActionPerformed(evt);
            }
        });
        jmDataset.add(miAddGene);

        miFuse.setText("Fuse gene entries");
        miFuse.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miFuseActionPerformed(evt);
            }
        });
        jmDataset.add(miFuse);

        miAlign.setText("Align strain sequences");
        miAlign.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miAlignActionPerformed(evt);
            }
        });
        jmDataset.add(miAlign);

        jmFix.setText("Fix");

        miRemoveInsertions.setText("Remove ORF disrupting insertions");
        miRemoveInsertions.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miRemoveInsertionsActionPerformed(evt);
            }
        });
        jmFix.add(miRemoveInsertions);

        miRemovePTC.setText("Remove premature terminal codons");
        miRemovePTC.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miRemovePTCActionPerformed(evt);
            }
        });
        jmFix.add(miRemovePTC);

        jmDataset.add(jmFix);
        jmDataset.add(jSeparator4);

        jmCT.setText("Codon table");
        jmCT.addMenuListener(new javax.swing.event.MenuListener() {
            public void menuCanceled(javax.swing.event.MenuEvent evt) {
            }
            public void menuDeselected(javax.swing.event.MenuEvent evt) {
            }
            public void menuSelected(javax.swing.event.MenuEvent evt) {
                jmCTMenuSelected(evt);
            }
        });
        jmDataset.add(jmCT);
        jmDataset.add(jSeparator7);

        miExit.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_X, java.awt.event.InputEvent.ALT_MASK));
        miExit.setText("Exit");
        miExit.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miExitActionPerformed(evt);
            }
        });
        jmDataset.add(miExit);

        jmMain.add(jmDataset);

        jmAnalyses.setText("Analyses");
        jmAnalyses.addMenuListener(new javax.swing.event.MenuListener() {
            public void menuCanceled(javax.swing.event.MenuEvent evt) {
            }
            public void menuDeselected(javax.swing.event.MenuEvent evt) {
            }
            public void menuSelected(javax.swing.event.MenuEvent evt) {
                jmAnalysesMenuSelected(evt);
            }
        });
        jmMain.add(jmAnalyses);

        jmFilters.setText("Filters");
        jmFilters.addMenuListener(new javax.swing.event.MenuListener() {
            public void menuCanceled(javax.swing.event.MenuEvent evt) {
            }
            public void menuDeselected(javax.swing.event.MenuEvent evt) {
            }
            public void menuSelected(javax.swing.event.MenuEvent evt) {
                jmFiltersMenuSelected(evt);
            }
        });
        jmMain.add(jmFilters);

        jmTools.setText("Tools");
        jmTools.addMenuListener(new javax.swing.event.MenuListener() {
            public void menuCanceled(javax.swing.event.MenuEvent evt) {
            }
            public void menuDeselected(javax.swing.event.MenuEvent evt) {
            }
            public void menuSelected(javax.swing.event.MenuEvent evt) {
                jmToolsMenuSelected(evt);
            }
        });

        jmSort.setText("Sort");

        miSortByName.setText("By gene name");
        miSortByName.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSortByNameActionPerformed(evt);
            }
        });
        jmSort.add(miSortByName);

        miSortByQuality.setText("By quality level");
        miSortByQuality.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSortByQualityActionPerformed(evt);
            }
        });
        jmSort.add(miSortByQuality);

        jmTools.add(jmSort);

        miQualityCheck.setText("Quality check");
        miQualityCheck.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miQualityCheckActionPerformed(evt);
            }
        });
        jmTools.add(miQualityCheck);

        miCreateCT.setText("Custom codon table");
        miCreateCT.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCreateCTActionPerformed(evt);
            }
        });
        jmTools.add(miCreateCT);

        miCoverage.setText("Coverage");
        miCoverage.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCoverageActionPerformed(evt);
            }
        });
        jmTools.add(miCoverage);

        miShowNames.setText("Show selected gene names");
        miShowNames.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miShowNamesActionPerformed(evt);
            }
        });
        jmTools.add(miShowNames);
        jmTools.add(jSeparator8);

        miScript.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.ALT_MASK));
        miScript.setText("Script");
        miScript.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miScriptActionPerformed(evt);
            }
        });
        jmTools.add(miScript);
        jmTools.add(jSeparator1);

        miSettings.setText("Settings");
        miSettings.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSettingsActionPerformed(evt);
            }
        });
        jmTools.add(miSettings);

        jmMain.add(jmTools);

        setJMenuBar(jmMain);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jSeparator3, javax.swing.GroupLayout.PREFERRED_SIZE, 6, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jViewsToolbar, javax.swing.GroupLayout.DEFAULT_SIZE, 434, Short.MAX_VALUE)
                    .addComponent(jpToolView, javax.swing.GroupLayout.DEFAULT_SIZE, 434, Short.MAX_VALUE))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jSeparator3, javax.swing.GroupLayout.DEFAULT_SIZE, 787, Short.MAX_VALUE)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jViewsToolbar, javax.swing.GroupLayout.PREFERRED_SIZE, 23, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jpToolView, javax.swing.GroupLayout.DEFAULT_SIZE, 731, Short.MAX_VALUE)
                .addContainerGap())
        );

        getAccessibleContext().setAccessibleName(null);

        pack();
    }// </editor-fold>//GEN-END:initComponents

private void miImportActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miImportActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            // Import the dataset.
            JFileChooser fc = new JFileChooser();
            String strLastOpenDir = sm.getSetting("", "LastOpenDir");
            if(strLastOpenDir!=null)
                fc.setCurrentDirectory(new File(strLastOpenDir));
            fc.setMultiSelectionEnabled(true);
            fc.setAcceptAllFileFilterUsed(false);
            for(PluginDetails pd:importers)
                fc.addChoosableFileFilter(pd.ed);
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            // Show OpenDialog.
            File[] files = null;
            if(fc.showOpenDialog(null)==JFileChooser.APPROVE_OPTION)
                files = fc.getSelectedFiles();
            // If no files are selected, do nothing.
            if(files==null)
                return;
            // Otherwise check whether the files exist.
            Vector<File> tmp = new Vector<File>();
            for(File f:files)
            {
                if(f.exists())
                    tmp.add(f);
            }
            // If no valid files are available, return.
            if(tmp.size()==0)
                return;
            // Update settings.
            String strLastDir = tmp.get(0).getParent();
            sm.addSetting("", "LastOpenDir", strLastDir);
            // Prepare an array of files and create the dataset.
            files = tmp.toArray(new File[1]);
            kernel.setSelectedIndices(getSelectedIndices());
            int index = ((ExtensionDetails)fc.getFileFilter()).getIndex();
            ErrorCode ec = kernel.importDataset(files, importers[index].strName, null);
            if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
            {
                String strErrDesc = kernel.getLastErrorString();
                String strErrMsg = null;
                if(strErrDesc==null || strErrDesc.isEmpty())
                    strErrMsg = String.format("<html>An unspecified error occured while importing the data with <b>%s</b></html>",
                                              importers[index].strName);
                else
                    strErrMsg = String.format("<html>An error occured while importing the data with <b>%s:</b><br>\t%s</html>",
                                              importers[index].strName, strErrDesc);
                JOptionPane.showMessageDialog(null, strErrMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            }
            updateGeneList();
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miImportActionPerformed

private void jmDatasetMenuSelected(javax.swing.event.MenuEvent evt) {//GEN-FIRST:event_jmDatasetMenuSelected
    boolean enabled = kernel.isDatasetLoaded();
    miExport.setEnabled(enabled);
    miUnload.setEnabled(enabled);
    jmFix.setEnabled(enabled);
    miFuse.setEnabled(enabled);
    miAlign.setEnabled(enabled);
}//GEN-LAST:event_jmDatasetMenuSelected

private void miExitActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miExitActionPerformed
    System.exit(0);
}//GEN-LAST:event_miExitActionPerformed

private void miUnloadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miUnloadActionPerformed
    kernel.unloadDataset();
    updateGeneList();
    for(int i=0;i<views.length;i++)
        views[i].displayGeneEntryInfo(null);
}//GEN-LAST:event_miUnloadActionPerformed

private void miExportActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miExportActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            // Code type.
            String strCodeType = sm.getSetting("", SettingsManager.CODETYPE);
            boolean bSimple = strCodeType.equalsIgnoreCase(SettingsManager.CT_SIMPLE);
            boolean bExtend = strCodeType.equalsIgnoreCase(SettingsManager.CT_EXTENDED);
            boolean bCompl  = strCodeType.equalsIgnoreCase(SettingsManager.CT_COMPLETE);
            // Save dialog.
            JFileChooser fc = new JFileChooser();
            String strLastSaveDir = sm.getSetting("", "LastSaveDir");
            if(strLastSaveDir!=null)
                fc.setCurrentDirectory(new File(strLastSaveDir));
            fc.setAcceptAllFileFilterUsed(false);
            for(PluginDetails pd:exporters)
            {
                if(bSimple)
                    fc.addChoosableFileFilter(pd.ed);
                else if(bExtend && pd.bSmd)
                    fc.addChoosableFileFilter(pd.ed);
                else if(bCompl && pd.bSad)
                    fc.addChoosableFileFilter(pd.ed);
            }
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            // Show OpenDialog.
            File file = null;
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            if(fc.showSaveDialog(null)==JFileChooser.APPROVE_OPTION)
                file = fc.getSelectedFile();
            // If no file is selected, do nothing.
            if(file==null)
                return;
            // If the returned file is a directory whereas a file is required,
            // add the file name <DATE_TIME>.
            if(file.isDirectory() && !fc.getFileFilter().getDescription().contains("(Folder only)"))
            {
                SimpleDateFormat sdf = new SimpleDateFormat("MMddyyyyHHmmss");
                Date date = new Date();
                file = new File(file.getAbsolutePath()+File.separator+sdf.format(date)+"."+
                                ((ExtensionDetails)(fc.getFileFilter())).getExtention()[0]);
            }
            // Save data set.
            sm.addSetting("", "LastSaveDir", file.getParent());
            int index = ((ExtensionDetails)fc.getFileFilter()).getIndex();
            ErrorCode ec = kernel.exportDataset(getSelectedIndices(), exporters[index].strName, file, null);
            if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
            {
                String strErrDesc = kernel.getLastErrorString();
                String strErrMsg = null;
                if(strErrDesc==null || strErrDesc.isEmpty())
                    strErrMsg = String.format("<html>An unspecified error occured while saving the data with <b>%s</b></html>",
                                              exporters[index].strName);
                else
                    strErrMsg = String.format("<html>An error occured while saving the data with <b>%s:</b><br>\t%s</html>",
                                              exporters[index].strName, strErrDesc);
                JOptionPane.showMessageDialog(null, strErrMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            }
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miExportActionPerformed

private void jmFiltersMenuSelected(javax.swing.event.MenuEvent evt) {//GEN-FIRST:event_jmFiltersMenuSelected
    boolean enabled = kernel.isDatasetLoaded();
    String strCodeType = sm.getSetting("", SettingsManager.CODETYPE);
    boolean bSimple = true;
    boolean bExtend = false;
    boolean bCompl = false;
    if(strCodeType!=null)
    {
        bSimple = strCodeType.equalsIgnoreCase(SettingsManager.CT_SIMPLE);
        bExtend = strCodeType.equalsIgnoreCase(SettingsManager.CT_EXTENDED);
        bCompl  = strCodeType.equalsIgnoreCase(SettingsManager.CT_COMPLETE);
    }
    int nFilters = filters.length;
    for(int i=0;i<nFilters;i++)
    {
        // If no dataset is loaded, disable ALL analyzers.
        if(!enabled)
            jmFilters.getItem(i).setEnabled(false);
        // Otherwise, if the current setting is complete code type, disable the analyzers which do
        // not support missing data and enable the others.
        else
        {
            // Enable the menu items depending on whether the corresponding
            // plugins support current code type.
            if(bSimple) // Enable all.
                jmFilters.getItem(i).setEnabled(true);
            else if(bExtend && filters[i].bSmd)
                jmFilters.getItem(i).setEnabled(true);
            else if(bCompl && filters[i].bSad )
                jmFilters.getItem(i).setEnabled(true);
            else
                jmFilters.getItem(i).setEnabled(false);
        }
    }
}//GEN-LAST:event_jmFiltersMenuSelected

private void jmiSelectActionPerformed1(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jmiSelectActionPerformed1
    DefaultTableModel model = (DefaultTableModel)lvGenes.getModel();
    for(int i=0;i<model.getRowCount();i++)
        model.setValueAt((Boolean)true, i, 0);
}//GEN-LAST:event_jmiSelectActionPerformed1

private void jmiInvertActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jmiInvertActionPerformed
    DefaultTableModel model = (DefaultTableModel)lvGenes.getModel();
    for(int i=0;i<model.getRowCount();i++)
        model.setValueAt(!(Boolean)model.getValueAt(i, 0), i, 0);
}//GEN-LAST:event_jmiInvertActionPerformed

private void lvGenesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_lvGenesMouseClicked
    if(evt.getClickCount()==1)
    {
        updateViews();
    }
    else if(evt.getClickCount()==2)
    {
        int i = lvGenes.getSelectedRow();
        if(i>-1)
            (new EditorForm()).editEntry(kernel.getDataset().getGeneEntry(i));
    }
}//GEN-LAST:event_lvGenesMouseClicked

private void jmiUnselectActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jmiUnselectActionPerformed
    DefaultTableModel model = (DefaultTableModel)lvGenes.getModel();
    for(int i=0;i<model.getRowCount();i++)
        model.setValueAt((Boolean)false, i, 0);
}//GEN-LAST:event_jmiUnselectActionPerformed

private void miScriptActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miScriptActionPerformed
    if(sf==null)
        sf = new ScriptForm(interpreter, kernel, (IWaitDialog)this.getGlassPane());
    sf.setVisible(true);
}//GEN-LAST:event_miScriptActionPerformed

private void miSettingsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSettingsActionPerformed
    (new SettingsDlg(kernel)).setVisible(true);
    kernel.updateSettings();
}//GEN-LAST:event_miSettingsActionPerformed

private void miCreateCTActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCreateCTActionPerformed
    if(ctc==null)
        ctc = new CodonTableCreator();
    CustomCodonTable ct = ctc.createCodonTable();
    if(ct!=null)
        kernel.addCodonTable(ct);
}//GEN-LAST:event_miCreateCTActionPerformed

private void jmAnalysesMenuSelected(javax.swing.event.MenuEvent evt) {//GEN-FIRST:event_jmAnalysesMenuSelected
    boolean enabled = kernel.isDatasetLoaded();
    String strCodeType = sm.getSetting("", SettingsManager.CODETYPE);
    boolean bSimple = true;
    boolean bExtend = false;
    boolean bCompl = false;
    if(strCodeType!=null)
    {
        bSimple = strCodeType.equalsIgnoreCase(SettingsManager.CT_SIMPLE);
        bExtend = strCodeType.equalsIgnoreCase(SettingsManager.CT_EXTENDED);
        bCompl  = strCodeType.equalsIgnoreCase(SettingsManager.CT_COMPLETE);
    }
    int nAnalyzers = analyzers.length;
    for(int i=0;i<nAnalyzers;i++)
    {
        // If no dataset is loaded, disable ALL analyzers.
        if(!enabled)
            jmAnalyses.getItem(i).setEnabled(false);
        // Otherwise, if the current setting is complete code type, disable the analyzers which do
        // not support missing data and enable the others.
        else
        {
            // Enable the menu items depending on whether the corresponding
            // plugins support current code type.
            if(bSimple) // Enable all.
                jmAnalyses.getItem(i).setEnabled(true);
            else if(bExtend && analyzers[i].bSmd)
                jmAnalyses.getItem(i).setEnabled(true);
            else if(bCompl && analyzers[i].bSad )
                jmAnalyses.getItem(i).setEnabled(true);
            else
                jmAnalyses.getItem(i).setEnabled(false);
        }
    }
}//GEN-LAST:event_jmAnalysesMenuSelected

private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
    Point p = this.getLocation();
    sm.addSetting("", WIN_X, Integer.toString(p.x));
    sm.addSetting("", WIN_Y, Integer.toString(p.y));
    Dimension d = this.getSize();
    sm.addSetting("", WIN_WIDTH, Integer.toString(d.width));
    sm.addSetting("", WIN_HEIGHT, Integer.toString(d.height));
    sm.addSetting("", WIN_STATE, Integer.toString(this.getExtendedState()));
    sm.addSetting("", SettingsManager.CODONTABLE, initData.ct.getName());
    sm.save("Settings.xml", false);
}//GEN-LAST:event_formWindowClosing

private void miCreateDpgpActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCreateDpgpActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            ErrorCode ec = kernel.importDPGPDataset(null);
            if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
            {
                String strErr = String.format("<html>An error occured while importing the DPGP data<br>\t%s</html>",
                        kernel.getLastErrorString());
                JOptionPane.showMessageDialog(null, strErr, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            }
            updateGeneList();
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miCreateDpgpActionPerformed

private void jmToolsMenuSelected(javax.swing.event.MenuEvent evt) {//GEN-FIRST:event_jmToolsMenuSelected
    boolean enabled = kernel.isDatasetLoaded();
    miSortByName.setEnabled(enabled);
    miSortByQuality.setEnabled(enabled);
    miQualityCheck.setEnabled(enabled);    
    miShowNames.setEnabled(enabled);
    miCoverage.setEnabled(enabled);
}//GEN-LAST:event_jmToolsMenuSelected

private void miSortByNameActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSortByNameActionPerformed
    kernel.sort();
    updateGeneList();
}//GEN-LAST:event_miSortByNameActionPerformed

private void miQualityCheckActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miQualityCheckActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            initData.wd.show(IWaitDialog.TYPE.Custom);
            kernel.runQualityCheck();
            updateGeneList();
            initData.wd.close();
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miQualityCheckActionPerformed

private void miRemoveInsertionsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miRemoveInsertionsActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            if(dsf==null)
                dsf = new DatasetFixer(initData);
            dsf.removeORFDisruptingInsertions(kernel.getDataset());
            initData.wd.show(IWaitDialog.TYPE.Kernel);
            kernel.runQualityCheck();
            updateGeneList();
            initData.wd.close();
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miRemoveInsertionsActionPerformed

private void miSortByQualityActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSortByQualityActionPerformed
    kernel.sortByQualityLevel();
    updateGeneList();
}//GEN-LAST:event_miSortByQualityActionPerformed

private void lvGenesKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_lvGenesKeyReleased
    int iCode = evt.getKeyCode();
    if(iCode==KeyEvent.VK_DOWN || iCode==KeyEvent.VK_UP)
        updateViews();
}//GEN-LAST:event_lvGenesKeyReleased

private void miRemovePTCActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miRemovePTCActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            if(dsf==null)
                dsf = new DatasetFixer(initData);
            dsf.maskPrematureTerminalCodons(kernel.getDataset());
            initData.wd.show(IWaitDialog.TYPE.Kernel);
            kernel.runQualityCheck();
            updateGeneList();
            initData.wd.close();
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miRemovePTCActionPerformed

private void miAlignActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miAlignActionPerformed
    IAligner[] aligners = kernel.listAligners();
    Dataset tmp = new Dataset();
    int[] indices = getSelectedIndices();
    if(indices.length==0)
        return;
    Dataset dataset = kernel.getDataset();
    for(int i:indices)
        tmp.addGene(dataset.getGeneEntry(i));
    AlignmentCreator ac = new AlignmentCreator(kernel, (IWaitDialog)this.getGlassPane());
    ac.showAlignmentDialog(aligners, tmp);    
}//GEN-LAST:event_miAlignActionPerformed

private void miFuseActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miFuseActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            ErrorCode ec = kernel.fuseGeneEntries(getSelectedIndices(), null);
            if(ec!=ErrorCode.Ok && ec!=ErrorCode.CancelledByUser)
            {
                String strErr = String.format("<html>An error occured while fusing the genes<br>\t%s</html>",
                        kernel.getLastErrorString());
                JOptionPane.showMessageDialog(null, strErr, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            }
            updateGeneList();
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miFuseActionPerformed

private void miAddGeneActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miAddGeneActionPerformed
    GeneEntry ge = (new AddGeneDialog()).createGeneEntry();
    if(ge==null)
        return;
    if(kernel.addGeneEntry(ge, false)==ErrorCode.ObjectExists)
    {
        String strMsg = String.format("<html>The dataset already contains a gene named <b>%s</b>. Add anyway?</html>");
        if(JOptionPane.showConfirmDialog(null, strMsg, Main.APPTITLE, JOptionPane.QUESTION_MESSAGE)==0)
        {
            kernel.addGeneEntry(ge, true);            
        }
    }
    updateGeneList();
}//GEN-LAST:event_miAddGeneActionPerformed

private void miShowNamesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miShowNamesActionPerformed
    int[] indices = getSelectedIndices();
    int nGenes = indices.length;
    if(nGenes==0)
    {
        JOptionPane.showMessageDialog(this, "No genes selected.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
        return;
    }
    Dataset ds = kernel.getDataset();
    String[] tmp = new String[nGenes];
    for(int i=0;i<nGenes;i++)
        tmp[i] = ds.getGeneEntry(indices[i]).getCommonName();
    (new ShowGenesDialog()).showNames(tmp);
}//GEN-LAST:event_miShowNamesActionPerformed

private void jmCTMenuSelected(javax.swing.event.MenuEvent evt) {//GEN-FIRST:event_jmCTMenuSelected
    String[] tmp = kernel.listCodonTables();
    ButtonGroup bg_ct = new ButtonGroup();
    jmCT.removeAll();
    for(String s:tmp)
    {
        JMenuItem mi = new JRadioButtonMenuItem(s);
        bg_ct.add(mi);
        mi.addActionListener(new CTMenuItemListener());
        jmCT.add(mi);
        if(s.equalsIgnoreCase(initData.ct.getName()))
            mi.setSelected(true);
    }
}//GEN-LAST:event_jmCTMenuSelected

private void miCoverageActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCoverageActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            IWaitDialog wd = kernel.getInitializationData().wd;
            wd.show(IWaitDialog.TYPE.Kernel);
            int[] indices = getSelectedIndices();
            int nGenes = indices.length;
            if(nGenes==0)
            {
                JOptionPane.showMessageDialog(null, "No genes selected.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
                return;
            }
            Dataset ds = kernel.getDataset();
            Dataset tmp = new Dataset();
            for(int i:indices)
                tmp.addGene(ds.getGeneEntry(i));
            Coverage cov = Coverage.calculateCoverage(tmp);
            wd.close();
            (new CoverageDetails()).displayCoverage(cov);
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_miCoverageActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator3;
    private javax.swing.JSeparator jSeparator4;
    private javax.swing.JSeparator jSeparator5;
    private javax.swing.JSeparator jSeparator6;
    private javax.swing.JSeparator jSeparator7;
    private javax.swing.JSeparator jSeparator8;
    private javax.swing.JToolBar jViewsToolbar;
    private javax.swing.JMenu jmAnalyses;
    private javax.swing.JMenu jmCT;
    private javax.swing.JMenu jmDataset;
    private javax.swing.JMenu jmFilters;
    private javax.swing.JMenu jmFix;
    private javax.swing.JMenu jmImport;
    private javax.swing.JMenuBar jmMain;
    private javax.swing.JMenu jmSort;
    private javax.swing.JMenu jmTools;
    private javax.swing.JMenuItem jmiInvert;
    private javax.swing.JMenuItem jmiSelect;
    private javax.swing.JMenuItem jmiUnselect;
    private javax.swing.JPanel jpToolView;
    private javax.swing.JPopupMenu jpmGenesMenu;
    private javax.swing.JTable lvGenes;
    private javax.swing.JMenuItem miAddGene;
    private javax.swing.JMenuItem miAlign;
    private javax.swing.JMenuItem miCoverage;
    private javax.swing.JMenuItem miCreateCT;
    private javax.swing.JMenuItem miCreateDpgp;
    private javax.swing.JMenuItem miExit;
    private javax.swing.JMenuItem miExport;
    private javax.swing.JMenuItem miFuse;
    private javax.swing.JMenuItem miImport;
    private javax.swing.JMenuItem miQualityCheck;
    private javax.swing.JMenuItem miRemoveInsertions;
    private javax.swing.JMenuItem miRemovePTC;
    private javax.swing.JMenuItem miScript;
    private javax.swing.JMenuItem miSettings;
    private javax.swing.JMenuItem miShowNames;
    private javax.swing.JMenuItem miSortByName;
    private javax.swing.JMenuItem miSortByQuality;
    private javax.swing.JMenuItem miUnload;
    // End of variables declaration//GEN-END:variables
}