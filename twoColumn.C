#include <iostream>
#include <string>
#include <fstream>

twoColumn(){
string root;
string plot;

cout <<" What is the .root file containing the histogram?"<<endl;
cin >> root;
TFile* file = new TFile(root.c_str());
file->ls();


while(1){
ofstream myfile;
string output;

cout <<" What is the name of the histogram to extract column data from?"<<endl;
cin.ignore(1000,'\n');
getline (cin,plot);
TH1D* h1 = (TH1D*)file->Get(plot.c_str());
double binWidth = h1->GetXaxis()->GetBinWidth(1);
cout <<" What is the name of the text file you want output (example.txt)?"<<endl;
cin >> output;
myfile.open(output.c_str());

cout<<endl<<endl<<endl;
myfile << "X-Axis Position (GeV)    Bin Content(dsigma/dx)" << endl;
int bins = h1->GetSize();
double x=0;
for (int i = 1;i<bins-1;i++){
x=binWidth+x;
myfile<<setw(12)<< x-binWidth/2<<"                    " << h1->GetBinContent(i)<< endl;
}

myfile <<endl<< " bins are " << h1->GetXaxis()->GetBinWidth(1) << " GeV wide "<<endl;
int entries =  h1->GetEntries();
myfile << " Number of entries = " << entries<< endl;
myfile << " Number of bins = " << bins-2<<endl;
myfile << " Note I am cutting off the first and last bins (over and under flow)"<<endl;
myfile.close();
cout << "Would you like to do another histogram from this file (y/n)?"<<endl;
string again;
cin >> again;
if(again == "n" || again == "N"){break;}
}
}
