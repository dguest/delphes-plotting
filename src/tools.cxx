#include "tools.hh"

#include "TFile.h"
#include "TKey.h"
#include "TCollection.h"

#include <stdexcept>
#include <iostream>
#include <string>
#include <regex>
#include <cassert>

using namespace std;

namespace TrackParam{
  enum trkParDef {D0=0, Z0, PHI, THETA, QOVERP};
  enum trkCovDef {D0D0=0, Z0D0, Z0Z0, PHID0, PHIZ0, PHIPHI, THETAD0, THETAZ0, THETAPHI, THETATHETA,
                  QOVERPD0, QOVERPZ0, QOVERPPHI, QOVERPTHETA, QOVERPQOVERP};

}


pair<int, int> pt_eta(const smatch& match) {
  return make_pair(stoi(match[1]), stoi(match[2]));
}

int dumpcov(string root_file) {
  TFile file_para(root_file.c_str(),"READ");
  if (!file_para.IsOpen() || file_para.IsZombie()) {
    throw std::runtime_error("bad file: " + string(root_file));
  }

  string regtail = "_ptbin([0-9]{2})_etabin([0-9]{2})";
  regex regex_cov_mat("^covmat" + regtail);
  regex regex_mean_vec("^meanvec" + regtail);

  TIter nextkey(file_para.GetListOfKeys());
  TKey *key;
  while ((key=(TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
    string name = key->GetName();
    std::smatch results;
    if (regex_search(name, results, regex_cov_mat)) {
      assert(results.size() == 3);
      auto bins = pt_eta(results);
      cout  << bins.first << " " << bins.second << endl;
      printf("at key:%s, object class:%s\n",key->GetName(),obj->ClassName());
    }
    delete obj;
  }
  return 0;
}

