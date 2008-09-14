// File created: 08/02/2008                          last modified: 09/05/2008
// Author: Christian Stratowa 08/02/2008

/*
 *******************************************************************************
 *
 * NOTE:
 * This macro is called by method "root.profilechart" of package "xps".
 * It allows to open a ROOT canvas to display a parallel coordinates plot
 *
 *******************************************************************************
 */

void macroDrawProfilePlot(const char *xps, const char *filename, const char *canvasname,
                          const char *treename, const char *varlist, const char *savename, 
                          Double_t min, Double_t max, Int_t aslog, Int_t globalscale,
                          Int_t candle, Int_t width, Int_t height)
{
   gSystem->Load("libGui");    //necessary for ROOT version >=5.15/09
//   gSystem->Load("libTreeViewer.so");  //necessary for TParallelCoord
   gSystem->Load(xps);

// create new manager
   XPlot *plotter = new XPlot("Plot");

// open root file
   plotter->Open(filename);

// add trees
   plotter->AddTree(treename);

   plotter->NewCanvas(canvasname, "title", 20, 20, width, height);

// draw in separate canvases
   plotter->DrawParallelCoord("", varlist, min, max, aslog, globalscale, candle);

// save figure and close canvas
   if (strcmp(savename, "") != 0) {
      plotter->SaveAs(savename);
      plotter->CloseCanvas("");
   }//if

// cleanup
   delete plotter;
}//macroDrawProfilePlot
