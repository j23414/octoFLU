#! /usr/bin/env python

""" Test octoFLU initial page
This is a test
To build the executable
pyinstaller -F octoFLU_gui.py
"""

import wx
import os
import os.path
#import subprocess
import sys
from shutil import which
from shutil import copyfile
import subprocess
from octoFLU import octoFLU

# Attempt at multi-threading, please someone help :(
#import time
#from threading import Thread
#from pubsub import pub

# === Globals
version="v1.0"

## === Threading so the GUI doesn't freeze
#"""
#website: http://www.blog.pythonlibrary.org/2010/05/22/wxpython-and-threads/
#website: https://pypubsub.readthedocs.io/en/v4.0.3/usage/usage_basic.html
#"""
#
#class TestThread(Thread):
#    """Test Worker Thread Class."""
# 
#    #----------------------------------------------------------------------
#    def __init__(self):
#        """Init Worker Thread Class."""
#        Thread.__init__(self)
#        self.start()    # start the thread
# 
#    #----------------------------------------------------------------------
#    def run(self):
#        """Run Worker Thread."""
#        # This is the code executing in the new thread.
#        for i in range(6):
#            time.sleep(10)
#            wx.CallAfter(self.postTime, i)
#        time.sleep(5)
#        wx.CallAfter(pub.sendMessage, "update", "Thread finished!")
# 
#    #----------------------------------------------------------------------
#    def postTime(self, amt):
#        """
#        Send time to GUI
#        """
#        amtOfTime = (amt + 1) * 10
#        pub.sendMessage("update", amtOfTime)

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
        self.makeblastdb_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = makeblastdb_path,
                                                message="Select a file",
                                                style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        blastn_path = self.config.Read("blastn")
        self.blastn_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = blastn_path,
                                           message="Select a file",
                                           style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        mafft_path = self.config.Read("mafft")
        self.mafft_fp = wx.FilePickerCtrl(pnl, wx.ID_ANY, path = mafft_path,
                                          message="Select a file",
                                           style = wx.FLP_OPEN|wx.FLP_USE_TEXTCTRL)

        fasttree_path = self.config.Read("fasttree")
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
        
#        ## === thread tutorial
#        self.displayLbl = wx.StaticText(pnl, label="Amount of time since thread started goes here")
#        self.btn = btn = wx.Button(pnl, label="Start Thread")
# 
#        btn.Bind(wx.EVT_BUTTON, self.onButton)
# 
#        sizer = wx.BoxSizer(wx.VERTICAL)
#        sizer.Add(self.displayLbl, 0, wx.ALL|wx.CENTER, 5)
#        sizer.Add(btn, 0, wx.ALL|wx.CENTER, 5)
#        #panel.SetSizer(sizer)
# 
#        # create a pubsub receiver
#        pub.subscribe(self.updateDisplay, "update")
#        ## === end thread tutorial
        
        fgs.AddMany([(makeblastdb_st),(self.makeblastdb_fp, wx.ID_ANY, wx.EXPAND),
                     (blastn_st), (self.blastn_fp, wx.ID_ANY, wx.EXPAND),
                     (mafft_st), (self.mafft_fp, wx.ID_ANY, wx.EXPAND),
                     (fasttree_st), (self.fasttree_fp, wx.ID_ANY, wx.EXPAND),
                     (reference_st), (self.reference_fp, wx.ID_ANY, wx.EXPAND),
                     (input_st), (self.input_fp, wx.ID_ANY, wx.EXPAND),
                     (helpmsg_st), (self.helpmsg_tc, wx.ID_ANY, wx.EXPAND),
                     (spacer_st),(self.run_bt, wx.ID_ANY,wx.EXPAND)]) #,
#                     (self.displayLbl),(self.btn, wx.ID_ANY, wx.EXPAND)])

        fgs.AddGrowableRow(6,1)
        fgs.AddGrowableCol(1,1)

        hbox.Add(fgs, proportion=1, flag=wx.ALL|wx.EXPAND, border=15)
        pnl.SetSizer(hbox)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run_bt)
        wx.Shell(command='echo "octoFLU shell initialized"')

        
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
        if(len(self.input_fp.GetPath())==0):
            self.helpmsg_tc.AppendText("  No input file provided... please select an input file.\n")
            return
#        outDir = self.input_fp.GetPath()+"_output"
#        
#        if not os.path.exists(outDir):
#            os.mkdir(outDir)
#        else:
#            print(outDir + " already exists")
            
        # Call octoFLU
        octoFLU(self.input_fp.GetPath(),
                BLASTN=self.blastn_fp.GetPath(),
                MAKEBLASTDB=self.makeblastdb_fp.GetPath(),
                MAFFT=self.mafft_fp.GetPath(),
                FASTTREE=self.fasttree_fp.GetPath(),
                reference=self.reference_fp.GetPath()
                )
        
#    def onButton(self, event):
#        """
#        Runs the thread
#        """
#        TestThread()
#        self.displayLbl.SetLabel("Thread started!")
#        btn = event.GetEventObject()
#        btn.Disable()
        
#    def updateDisplay(self, msg):
#        """
#        Receives data from thread and updates the display
#        """
#        t = msg.data
#        if isinstance(t, int):
#            self.displayLbl.SetLabel("Time since thread started: %s seconds" % t)
#        else:
#            self.displayLbl.SetLabel("%s" % t)
#            self.btn.Enable()
            
            
if __name__ == '__main__':
    app = wx.App()
    frm = MainFrame(None, title="octoFLU " + version, size=(600,500))
    frm.Show()
    app.MainLoop()
    
