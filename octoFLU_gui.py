#! /usr/bin/env python

""" Test octoFLU initial page
This is a test
"""

import wx
import os
import os.path
import subprocess
import sys
from shutil import which
from shutil import copyfile
import subprocess
from octoFLU import octoFLU

# === Globals
version="v1.0"

class MainFrame(wx.Frame):
    """
    This is a basic frame
    """

    def __init__(self,*args,**kw):
        super(MainFrame, self).__init__(*args, **kw)
        self.config = wx.Config("octoFLU")
        self.InitUI()
        self.Centre()

        
    def InitUI(self):
        # Increase font size
        font = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)
        font.SetPointSize(9)

        # Panel and sizers
        pnl = wx.Panel(self)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        fgs = wx.FlexGridSizer(cols=2, vgap = 9, hgap= 25)

        # text Controls
        makeblastdb_st = wx.StaticText(pnl, label="makeblastdb path:")
        blastn_st = wx.StaticText(pnl, label="blastn path:")
        mafft_st = wx.StaticText(pnl, label="MAFFT path:")
        fasttree_st = wx.StaticText(pnl, label="FastTree path:")
        reference_st = wx.StaticText(pnl, label="Reference file:")
        input_st = wx.StaticText(pnl, label="Input file:")
        helpmsg_st = wx.StaticText(pnl, label="Help messages:")
        spacer_st = wx.StaticText(pnl, label="")

        makeblastdb_st.SetFont(font)
        blastn_st.SetFont(font)
        mafft_st.SetFont(font)
        fasttree_st.SetFont(font)
        reference_st.SetFont(font)
        input_st.SetFont(font)
        helpmsg_st.SetFont(font)

        makeblastdb_path = self.config.Read("makeblastdb")
#        if(len(makeblastdb_path)<1):
#            makeblastdb_path = os.popen("which makeblastdb").read()
        self.makeblastdb_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = makeblastdb_path,
                                                message="Select a file",
                                                style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        blastn_path = self.config.Read("blastn")
#        if(len(blastn_path)<1):
#            blastn_path = os.popen("which blastn").read()
        self.blastn_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = blastn_path,
                                           message="Select a file",
                                           style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        mafft_path = self.config.Read("mafft")
#        if(len(mafft_path)<1):
#            mafft_path = os.popen("which mafft").read()
        self.mafft_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = mafft_path,
                                          message="Select a file",
                                           style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        fasttree_path = self.config.Read("fasttree")
#        if(len(fasttree_path)<1):
#            fasttree_path = os.popen("which fasttree").read()
        self.fasttree_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = fasttree_path,
                                             message="Select a file",
                                           style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        self.reference_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = self.config.Read("reference"),
                                              message="Select a file",
                                           style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)
        
        self.input_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = wx.EmptyString,
                                          message="Select a file",
                                          style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)
        self.helpmsg_tc = wx.TextCtrl(pnl, style=wx.TE_MULTILINE)
        self.helpmsg_tc.AppendText("octoFLU initialized!\nPlease link \n\t * your executables \n\t * the octoFLU provided reference fasta file \n\t * and your query fasta file above:\n")
        
        sys.stdout=self.helpmsg_tc
#        sys.stderr=self.helpmsg_tc

        self.run_bt = wx.Button(pnl, wx.ID_ANY, label="Classify influenza genes")
        
        fgs.AddMany([(makeblastdb_st),(self.makeblastdb_fp, wx.ID_ANY, wx.EXPAND),
                     (blastn_st), (self.blastn_fp, wx.ID_ANY, wx.EXPAND),
                     (mafft_st), (self.mafft_fp, wx.ID_ANY, wx.EXPAND),
                     (fasttree_st), (self.fasttree_fp, wx.ID_ANY, wx.EXPAND),
                     (reference_st), (self.reference_fp, wx.ID_ANY, wx.EXPAND),
                     (input_st), (self.input_fp, wx.ID_ANY, wx.EXPAND),
                     (helpmsg_st), (self.helpmsg_tc, wx.ID_ANY, wx.EXPAND),
                     (spacer_st),(self.run_bt, wx.ID_ANY,wx.EXPAND)])

        fgs.AddGrowableRow(6,1)
        fgs.AddGrowableCol(1,1)

        hbox.Add(fgs, proportion=1, flag=wx.ALL|wx.EXPAND, border=15)
        pnl.SetSizer(hbox)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run_bt)
        wx.Shell(command='echo "Hello world"')



        
    def onFilePicker(self, event):
        self.resetOnOpen(event)
        path = self.filePicker.GetPath()
        if not os.path.isfile(path):
            return
        self.openFile(event, path)
        self.modifyHistory(event,path)

    def OnClose(self, event):
        # Save program paths in config file
        self.config.Write("makeblastdb", self.makeblastdb_fp.GetPath())
        self.config.Write("blastn", self.blastn_fp.GetPath())
        self.config.Write("mafft", self.mafft_fp.GetPath())
        self.config.Write("fasttree", self.fasttree_fp.GetPath())
        self.config.Write("reference", self.reference_fp.GetPath())
        # Exit program
        event.Skip()

    def OnRun(self, event):
        outDir = self.input_fp.GetPath()+"_output"
        
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        else:
            print(outDir + " already exists")
            
        # Call octoFLU
        octoFLU(self.input_fp.GetPath(),
                BLASTN=self.blastn_fp.GetPath(),
                MAKEBLASTDB=self.makeblastdb_fp.GetPath(),
                MAFFT=self.mafft_fp.GetPath(),
                FASTTREE=self.fasttree_fp.GetPath(),
                reference=self.reference_fp.GetPath()
                )

#        wx.Shell(command="echo running...")
#        wx.Shell(command="echo " + self.mafft_fp.GetPath()
        # ==== Create your Blast Database
#        subprocess.call([self.makeblastdb_fp.GetPath(),"-in",self.reference_fp.GetPath(),"-dbtype","nucl"])
        # ==== Search your blast database
#        subprocess.call([self.blastn_fp.GetPath(),"-db",self.reference_fp.GetPath(),"-query",self.input_fp.GetPath(),"-num_alignments","1","-outfmt","6", "-out",outDir+"/blast_output.txt"])
#        wx.Shell(command= "\"" + self.mafft_fp.GetPath() + "\"")
        

if __name__ == '__main__':
    app = wx.App()
    frm = MainFrame(None, title="octoFLU " + version, size=(600,500))
    frm.Show()
    app.MainLoop()
    
