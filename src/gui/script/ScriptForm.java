/*
    File:
        ScriptForm.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Script Editor window.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.script;

import algorithms.Time;
import dpgp.DPGPImporter;
import gui.IWaitDialog;
import java.awt.Cursor;
import kernel.PluginDetails;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;
import kernel.ErrorCode;
import kernel.Kernel;
import kernel.Main;
import kernel.SettingsManager;
import kernel.script.ScriptInterpreter;
import plugin.PluginType;


public class ScriptForm extends javax.swing.JFrame 
{
    private final String FILE_DESC = "GeneAnalyzer Script (*.gas)";
    private final String FILE_EXT  = "gas";
    
    private String strFileName              = null;
    private ScriptInterpreter interpreter   = null;
    private Kernel kernel                   = null;
    private SyntaxHighlighter shl           = null;
    private IWaitDialog wd                  = null;
    private SettingsManager sm              = null;

    
    /**
     *  Constructs the ScriptForm and the corresponding ScriptInterpreter.
     * 
     *  @param kernel
     */
    public ScriptForm(ScriptInterpreter interpreter, Kernel kernel, IWaitDialog wd)
    {
        super("Script Editor");
        initComponents();
        setLocationRelativeTo(null);
        this.interpreter = interpreter;
        this.kernel = kernel;
        shl = new SyntaxHighlighter(jtpCode);
        new Thread(shl).start();
        this.wd = wd;
        this.sm = kernel.getInitializationData().sm;
    }    
    
    /**
     *  Saves the current script to file.
     * 
     *  @param strFilename
     */
    private void saveToFile(String strFilename)
    {
        // If the filename is not specified, show SaveFileDialog.
        if(strFilename==null)
        {
            setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            JFileChooser fc = new JFileChooser();
            String strLastScriptDir = sm.getSetting("", "LastScriptSaveDir");
            if(strLastScriptDir!=null)
                fc.setCurrentDirectory(new File(strLastScriptDir));
            fc.addChoosableFileFilter(
                    new FileNameExtensionFilter(FILE_DESC, FILE_EXT));
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
            {
                strFilename = fc.getSelectedFile().getAbsolutePath();
                String strLastDir = fc.getSelectedFile().getParent();
                sm.addSetting("", "LastScriptSaveDir", strLastDir);
            }
            else
                return;
        }
        // Update global filename.
        this.strFileName = strFilename;
        // Save the file.
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter(new File(strFileName)));
            out.print(jtpCode.getText());
            out.close();
        }
        catch (IOException e)
        {
            JOptionPane.showMessageDialog(null, 
                                          "An I/O error occured while saving the script", 
                                          Main.APPTITLE, 
                                          JOptionPane.ERROR_MESSAGE);
            return;
        }
    }
    
    /**
     *  Adds the intruction to the code.
     * 
     *  @param strInstruction
     */
    private void addInstruction(String strInstruction)
    {
        if(jtpCode.getText().isEmpty())
            jtpCode.setText(strInstruction+"\n");
        else
        {
            if(jtpCode.getText().endsWith("\n")) 
                jtpCode.setText(jtpCode.getText()+strInstruction+"\n");
            else
                jtpCode.setText(jtpCode.getText()+"\n"+strInstruction+"\n");
        }
    }
    
    
    /**
     *  Shows the dialog with all available plugins of the specified type
     *  and allows to select one of them. If the user clicks CANCEL in the
     *  selection dialog, null is returned. 
     * 
     *  @param type
     *  @return
     */
    private PluginDetails getPluginDetails(PluginType type)
    {
        PluginSelectionForm psf = new PluginSelectionForm();
        PluginDetails[] pd = kernel.getPluginsDetails(type);
        if(pd.length==0)
        {
            JOptionPane.showMessageDialog(null, 
                                  "No appropriate plugins for this instruction found", 
                                  Main.APPTITLE, 
                                  JOptionPane.ERROR_MESSAGE);
            return null;
        }
        int index = psf.getPluginIndex(pd);
        if(index==-1)
            return null;
        else
            return pd[index];
    }

    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPopupMenu1 = new javax.swing.JPopupMenu();
        miLoad = new javax.swing.JMenuItem();
        miDpgpLoad = new javax.swing.JMenuItem();
        miAnalyze = new javax.swing.JMenuItem();
        miFilter = new javax.swing.JMenuItem();
        miUnload = new javax.swing.JMenuItem();
        miSave = new javax.swing.JMenuItem();
        miCodonTable = new javax.swing.JMenuItem();
        miSelectAll = new javax.swing.JMenuItem();
        miSortByName = new javax.swing.JMenuItem();
        miSortByQuality = new javax.swing.JMenuItem();
        miInvertSelection = new javax.swing.JMenuItem();
        miExecute = new javax.swing.JMenuItem();
        jSeparator3 = new javax.swing.JSeparator();
        miComm = new javax.swing.JMenuItem();
        jToolBar1 = new javax.swing.JToolBar();
        jbNew = new javax.swing.JButton();
        jbClear = new javax.swing.JButton();
        jbOpen = new javax.swing.JButton();
        jbSave = new javax.swing.JButton();
        jbSaveas = new javax.swing.JButton();
        jSeparator2 = new javax.swing.JToolBar.Separator();
        jbAdd = new javax.swing.JButton();
        jSeparator1 = new javax.swing.JToolBar.Separator();
        jbRun = new javax.swing.JButton();
        jScrollPane1 = new javax.swing.JScrollPane();
        jtpCode = new javax.swing.JTextPane();
        jLabel1 = new javax.swing.JLabel();
        jlRunTime = new javax.swing.JLabel();

        miLoad.setText("LOAD");
        miLoad.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miLoadActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miLoad);

        miDpgpLoad.setText("LOADDPGP");
        miDpgpLoad.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miDpgpLoadActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miDpgpLoad);

        miAnalyze.setText("ANALYZE");
        miAnalyze.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miAnalyzeActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miAnalyze);

        miFilter.setText("FILTER");
        miFilter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miFilterActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miFilter);

        miUnload.setText("UNLOAD");
        miUnload.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miUnloadActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miUnload);

        miSave.setText("SAVE");
        miSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSaveActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miSave);

        miCodonTable.setText("CODONTABLE");
        miCodonTable.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCodonTableActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miCodonTable);

        miSelectAll.setText("SELECT ALL");
        miSelectAll.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSelectAllActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miSelectAll);

        miSortByName.setText("SORT BY NAME");
        miSortByName.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSortByNameActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miSortByName);

        miSortByQuality.setText("SORT BY QUALITY");
        miSortByQuality.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miSortByQualityActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miSortByQuality);

        miInvertSelection.setText("INVERT SELECTION");
        miInvertSelection.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miInvertSelectionActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miInvertSelection);

        miExecute.setText("EXECUTE");
        miExecute.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miExecuteActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miExecute);
        jPopupMenu1.add(jSeparator3);

        miComm.setText("Comment");
        miComm.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCommActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miComm);

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Script Editor");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosed(java.awt.event.WindowEvent evt) {
                formWindowClosed(evt);
            }
        });

        jToolBar1.setFloatable(false);
        jToolBar1.setRollover(true);

        jbNew.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_new.png"))); // NOI18N
        jbNew.setText("New");
        jbNew.setToolTipText("Create a new script");
        jbNew.setFocusable(false);
        jbNew.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbNew.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbNew.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbNewActionPerformed(evt);
            }
        });
        jToolBar1.add(jbNew);

        jbClear.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_clear.png"))); // NOI18N
        jbClear.setText("Clear");
        jbClear.setToolTipText("Clear the code window");
        jbClear.setFocusable(false);
        jbClear.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbClear.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbClear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbClearActionPerformed(evt);
            }
        });
        jToolBar1.add(jbClear);

        jbOpen.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_open.png"))); // NOI18N
        jbOpen.setText("Open");
        jbOpen.setFocusable(false);
        jbOpen.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbOpen.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbOpen.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbOpenActionPerformed(evt);
            }
        });
        jToolBar1.add(jbOpen);

        jbSave.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_save.png"))); // NOI18N
        jbSave.setText("Save");
        jbSave.setEnabled(false);
        jbSave.setFocusable(false);
        jbSave.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbSave.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbSaveActionPerformed(evt);
            }
        });
        jToolBar1.add(jbSave);

        jbSaveas.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_saveas.png"))); // NOI18N
        jbSaveas.setText("Save as...");
        jbSaveas.setEnabled(false);
        jbSaveas.setFocusable(false);
        jbSaveas.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbSaveas.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbSaveas.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbSaveasActionPerformed(evt);
            }
        });
        jToolBar1.add(jbSaveas);
        jToolBar1.add(jSeparator2);

        jbAdd.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_add.png"))); // NOI18N
        jbAdd.setText("Add");
        jbAdd.setToolTipText("Add an instruction backbone");
        jbAdd.setFocusable(false);
        jbAdd.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbAdd.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbAdd.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAddActionPerformed(evt);
            }
        });
        jToolBar1.add(jbAdd);
        jToolBar1.add(jSeparator1);

        jbRun.setIcon(new javax.swing.ImageIcon(getClass().getResource("/img/ga_run.png"))); // NOI18N
        jbRun.setText("Run");
        jbRun.setToolTipText("Run the script");
        jbRun.setEnabled(false);
        jbRun.setFocusable(false);
        jbRun.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbRun.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbRun.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRunActionPerformed(evt);
            }
        });
        jToolBar1.add(jbRun);

        jtpCode.addCaretListener(new javax.swing.event.CaretListener() {
            public void caretUpdate(javax.swing.event.CaretEvent evt) {
                jtpCodeCaretUpdate(evt);
            }
        });
        jScrollPane1.setViewportView(jtpCode);

        jLabel1.setText("Last run:");

        jlRunTime.setText("     ");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 482, Short.MAX_VALUE)
                    .addComponent(jToolBar1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addGap(18, 18, 18)
                        .addComponent(jlRunTime)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jToolBar1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 313, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jlRunTime))
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

private void jbClearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbClearActionPerformed
    jtpCode.setText("");
}//GEN-LAST:event_jbClearActionPerformed

private void jbSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSaveActionPerformed
    saveToFile(strFileName);
}//GEN-LAST:event_jbSaveActionPerformed

private void jbSaveasActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbSaveasActionPerformed
    saveToFile(null);
}//GEN-LAST:event_jbSaveasActionPerformed

private void jbNewActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbNewActionPerformed
    if(!jtpCode.getText().isEmpty())
    {
        if(!jtpCode.getText().isEmpty())
        {
            if(JOptionPane.showConfirmDialog(null, 
                                             "Save changes?", 
                                             Main.APPTITLE, 
                                             JOptionPane.YES_NO_OPTION, 
                                             JOptionPane.QUESTION_MESSAGE)==0)
                saveToFile(strFileName);
        }
    }
    jtpCode.setText("");
    strFileName = null;
}//GEN-LAST:event_jbNewActionPerformed

private void jbRunActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRunActionPerformed
    Runnable r = new Runnable()
    {
        public void run()
        {
            long d = System.currentTimeMillis();
            wd.show(IWaitDialog.TYPE.Kernel);
            setVisible(false);
            if(interpreter.runScript(jtpCode.getText())!=ErrorCode.Ok)
            {
                JOptionPane.showMessageDialog(null, interpreter.getLastErrorString(),
                        Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            }
            wd.close();
            setVisible(true);
            d = System.currentTimeMillis() - d;
            Time time = Time.constructFromMillis(d);
            jlRunTime.setText(time.toString());
        }
    };
    new Thread(r).start();
}//GEN-LAST:event_jbRunActionPerformed

private void jbOpenActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbOpenActionPerformed
    if(!jtpCode.getText().isEmpty())
    {
        if(JOptionPane.showConfirmDialog(null, 
                                             "Save changes?", 
                                             Main.APPTITLE, 
                                             JOptionPane.YES_NO_OPTION, 
                                             JOptionPane.QUESTION_MESSAGE)==0)
            saveToFile(strFileName);
    }
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    // Show selection dialog.
    JFileChooser fc = new JFileChooser();
    String strLastScriptDir = sm.getSetting("", "LastScriptOpenDir");
    if(strLastScriptDir!=null)
        fc.setCurrentDirectory(new File(strLastScriptDir));
    fc.setAcceptAllFileFilterUsed(false);
    fc.addChoosableFileFilter(new FileNameExtensionFilter(FILE_DESC, FILE_EXT));
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
    if(fc.showOpenDialog(null)==JFileChooser.APPROVE_OPTION)
        strFileName = fc.getSelectedFile().getAbsolutePath();
    else
        return;
    // Load the content.
    try
    {
        BufferedReader in = new BufferedReader(new InputStreamReader(
            new FileInputStream(new File(strFileName))));
        StringBuffer sb = new StringBuffer();
        String strLine;
        while((strLine = in.readLine())!=null)
        {
            sb.append(strLine);
            sb.append("\n");            
        }
        jtpCode.setText(sb.toString());
        // Save the directory setting.
        String strLastDir = fc.getSelectedFile().getParent();
        sm.addSetting("", "LastScriptOpenDir", strLastDir);
    }
    catch(Exception e)
    {
        JOptionPane.showMessageDialog(null, 
                                      "An I/O error occured while reading the file",
                                      Main.APPTITLE, 
                                      JOptionPane.ERROR_MESSAGE);
        strFileName = null;
        return;
    }
}//GEN-LAST:event_jbOpenActionPerformed

private void jbAddActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbAddActionPerformed
    jPopupMenu1.show(this, jbAdd.getX(), jToolBar1.getY()+jToolBar1.getHeight());
}//GEN-LAST:event_jbAddActionPerformed

private void miUnloadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miUnloadActionPerformed
    addInstruction("UNLOAD;");
}//GEN-LAST:event_miUnloadActionPerformed

private void miLoadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miLoadActionPerformed
    PluginDetails pd = getPluginDetails(PluginType.IMPORTER);
    if(pd!=null)
    {
        addInstruction(String.format("LOAD name=\"%s\" source=\"<FILENAME>\" params=\"%s\";", 
                pd.strName, pd.strParamStr));
    }
}//GEN-LAST:event_miLoadActionPerformed

private void miAnalyzeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miAnalyzeActionPerformed
    PluginDetails pd = getPluginDetails(PluginType.ANALYZER);
    if(pd!=null)
    {
        addInstruction(String.format("ANALYZE name=\"%s\" params=\"%s\";", 
                pd.strName, pd.strParamStr));
    }
}//GEN-LAST:event_miAnalyzeActionPerformed

private void miFilterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miFilterActionPerformed
    PluginDetails pd = getPluginDetails(PluginType.FILTER);
    if(pd!=null)
    {
        addInstruction(String.format("FILTER name=\"%s\" params=\"%s\";", 
                pd.strName, pd.strParamStr));
    }
}//GEN-LAST:event_miFilterActionPerformed

private void miSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSaveActionPerformed
    PluginDetails pd = getPluginDetails(PluginType.EXPORTER);
    if(pd!=null)
    {
        addInstruction(String.format("SAVE name=\"%s\" dest=\"<FILENAME>\" params=\"%s\";", 
                pd.strName, pd.strParamStr));
    }
}//GEN-LAST:event_miSaveActionPerformed

private void miCommActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCommActionPerformed
    addInstruction("// Comment text here.");
}//GEN-LAST:event_miCommActionPerformed

private void miDpgpLoadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miDpgpLoadActionPerformed
    addInstruction(String.format("DPGPLOAD params=\"%s\";", DPGPImporter.getParamString()));
}//GEN-LAST:event_miDpgpLoadActionPerformed

private void miCodonTableActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCodonTableActionPerformed
    String strName = (new CodonTableSelectionForm()).getCodonTableName(kernel.listCodonTables());
    if(strName!=null)
        addInstruction(String.format("SYSTEM CODONTABLE \"%s\";", strName));
}//GEN-LAST:event_miCodonTableActionPerformed

private void jtpCodeCaretUpdate(javax.swing.event.CaretEvent evt) {//GEN-FIRST:event_jtpCodeCaretUpdate
    if(jtpCode.getText().isEmpty())
    {
        jbSave.setEnabled(false);
        jbRun.setEnabled(false);
        jbSaveas.setEnabled(false);
    }
    else
    {
        jbSave.setEnabled(true);
        jbRun.setEnabled(true);
        jbSaveas.setEnabled(true);
    }
}//GEN-LAST:event_jtpCodeCaretUpdate

private void formWindowClosed(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosed

}//GEN-LAST:event_formWindowClosed

private void miSelectAllActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSelectAllActionPerformed
    addInstruction("SYSTEM SELECTALL;");
}//GEN-LAST:event_miSelectAllActionPerformed

private void miSortByNameActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSortByNameActionPerformed
    addInstruction("SYSTEM SORT \"BY NAME\";");
}//GEN-LAST:event_miSortByNameActionPerformed

private void miSortByQualityActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miSortByQualityActionPerformed
    addInstruction("SYSTEM SORT \"BY QUALITY\";");
}//GEN-LAST:event_miSortByQualityActionPerformed

private void miInvertSelectionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miInvertSelectionActionPerformed
    addInstruction("SYSTEM INVERTSELECTION;");
}//GEN-LAST:event_miInvertSelectionActionPerformed

private void miExecuteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miExecuteActionPerformed
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    JFileChooser fc = new JFileChooser();
    String strLastExecDir = sm.getSetting("", "LastExecuteOpenDir");
    if(strLastExecDir!=null)
        fc.setCurrentDirectory(new File(strLastExecDir));
    fc.setAcceptAllFileFilterUsed(true);
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
    String strCmd = null;
    if(fc.showOpenDialog(null)==JFileChooser.APPROVE_OPTION)
        strCmd = fc.getSelectedFile().getAbsolutePath();
    else
        return;
    strLastExecDir = fc.getSelectedFile().getParent();
    sm.addSetting("", "LastExecuteOpenDir", strLastExecDir);
    addInstruction(String.format("SYSTEM EXECUTE NOWAIT \"%s\";", strCmd));
}//GEN-LAST:event_miExecuteActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel1;
    private javax.swing.JPopupMenu jPopupMenu1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JToolBar.Separator jSeparator1;
    private javax.swing.JToolBar.Separator jSeparator2;
    private javax.swing.JSeparator jSeparator3;
    private javax.swing.JToolBar jToolBar1;
    private javax.swing.JButton jbAdd;
    private javax.swing.JButton jbClear;
    private javax.swing.JButton jbNew;
    private javax.swing.JButton jbOpen;
    private javax.swing.JButton jbRun;
    private javax.swing.JButton jbSave;
    private javax.swing.JButton jbSaveas;
    private javax.swing.JLabel jlRunTime;
    private javax.swing.JTextPane jtpCode;
    private javax.swing.JMenuItem miAnalyze;
    private javax.swing.JMenuItem miCodonTable;
    private javax.swing.JMenuItem miComm;
    private javax.swing.JMenuItem miDpgpLoad;
    private javax.swing.JMenuItem miExecute;
    private javax.swing.JMenuItem miFilter;
    private javax.swing.JMenuItem miInvertSelection;
    private javax.swing.JMenuItem miLoad;
    private javax.swing.JMenuItem miSave;
    private javax.swing.JMenuItem miSelectAll;
    private javax.swing.JMenuItem miSortByName;
    private javax.swing.JMenuItem miSortByQuality;
    private javax.swing.JMenuItem miUnload;
    // End of variables declaration//GEN-END:variables
}
