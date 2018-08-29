#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "TLorentzVector.h"


void PRINT4VECTOR( TLorentzVector v, bool doPxPyPzE ){
	

	cout << " --------------------- " << endl;
	if( doPxPyPzE ){
		cout << " Px = " << v.Px() << endl;
		cout << " Py = " << v.Py() << endl;
		cout << " Pz = " << v.Pz() << endl;
		cout << " Mass = " << v.M() << endl;
		cout << " E = " << v.E() << endl;
	}


}