#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
using namespace std;
using namespace erhic;

void checkFermi(){

	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_devK/k-new_v1.root" ); // Wild cards are allowed e.g. tree.Add("*.root" );

	EventPythia* event(NULL);// = new EventPythia;

	// EventBase* event(NULL);
	// EventBeagle* event_beagle(NULL);

	tree->SetBranchAddress("event", &event ); // Note &event, not event.

	TH1D* hist_pxf = new TH1D("hist_pxf","",1000,0,1);
	TH1D* hist_pyf = new TH1D("hist_pyf","",1000,0,1);
	TH1D* hist_pzf = new TH1D("hist_pzf","",1000,0,1);

	tree->Draw("pxf>>hist_pxf");
	tree->Draw("pxf>>hist_pyf");
	tree->Draw("pxf>>hist_pzf");

	hist_pxf->Draw();

}