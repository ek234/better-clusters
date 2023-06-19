/*
 *   This file generates ground truth jet clusters.
 */

#include "./imports.h"

void run_preamble(int argc, char *argv[], ofstream &fout) {
  ClusterSequence::print_banner();
  cout << endl;
  gSystem->Load("libDelphes");

  if (argc < 2) {
    cerr << "error. usage: " << argv[0] << " output_file" << endl;
    exit(EXIT_FAILURE);
  }

  fout.open(argv[1], fstream::out);
  while (!fout.is_open())
    ;
  fout << "event,"
       // << "parent_pid,"
       << "mass,"
       // << "charge,"
       << "pt,"
       << "energy,"
       << "px,"
       << "py,"
       << "pz,"
       << "girth,"
       << "variance,"
       << "skew,"
       << "kurt,"
       << "var_fake,"
       << "skew_fake,"
       << "kurt_fake,"
       << "tau1,"
       << "tau2,"
       << "tau3,"
       << "tau4,"
       << "tau5,"
       << "tau21,"
       << "tau32,"
       << "tau43,"
       << "tau54,"
       << "angularity_2,"
       << "angularity_3,"
       << "angularity_4,"
       << "angularity_5,"
       << "angularity_9" << endl;
  fout.flush();
}

double correctedDiff(double phi_a, double phi_b) {
  double diff = abs(phi_a - phi_b);
  return min(diff, 2 * 3.14159265 - diff);
}
double deltaR_pj(PseudoJet &a, PseudoJet &b) {
  return sqrt(pow(a.eta() - b.eta(), 2.) +
              pow(correctedDiff(a.phi(), b.phi()), 2.));
}

double calculateGirth(PseudoJet &jet) {
  double girth = 0;
  for (auto &element : jet.constituents())
    girth += element.pt() * deltaR_pj(element, jet);
  girth /= jet.pt();
  return girth;
}

vector<double> calculateHigherMoments(PseudoJet &jet, double girth,
                                      bool mistake = false) {
  double deviation_2 = 0;
  double deviation_3 = 0;
  double deviation_4 = 0;

  for (auto &element : jet.constituents()) {
    double diff;
    if (not mistake)
      diff = element.pt() * (deltaR_pj(element, jet) - girth);
    else
      diff = (element.pt() * deltaR_pj(element, jet)) - girth;

    double diff2 = diff * diff;
    double diff3 = diff2 * diff;
    double diff4 = diff2 * diff2;

    deviation_2 += diff2;
    deviation_3 += diff2;
    deviation_4 += diff4;
  }

  double s = deviation_2 / jet.pt();
  double skew = (deviation_3 / jet.pt()) / (powf64(s, 1.5));
  double kurt = (deviation_4 / jet.pt()) / (powf64(s, 2));

  return {s, skew, kurt};
}

double calculateAngularity(PseudoJet &jet, double a) {
  double angularity = 0;
  for (auto &element : jet.constituents())
    angularity += element.pt() * exp((a - 1) * abs(element.rap() - jet.rap()));
  angularity /= 2 * jet.E();
  return angularity;
}

vector<double> calculateNsub(PseudoJet &pjet, double beta = 1.0) {
  // beta = 1:  aka broadening/girth/width measure
  //   wta_kt_axes are approximately the same as minimizing beta = 1 measure
  //
  // beta = 2:  aka thrust/mass measure
  //   kt_axes are approximately the same as minimizing beta = 2 measure
  //
  // N.B. The minimization routines are only valid for 1 < beta < 3.

  UnnormalizedMeasure measureSpec(beta);
  OnePass_WTA_KT_Axes axisMode; // Why not AntiKT

  Nsubjettiness nSub1(1, axisMode, measureSpec);
  Nsubjettiness nSub2(2, axisMode, measureSpec);
  Nsubjettiness nSub3(3, axisMode, measureSpec);
  Nsubjettiness nSub4(4, axisMode, measureSpec);
  Nsubjettiness nSub5(5, axisMode, measureSpec);
  NsubjettinessRatio nSub21(2, 1, axisMode, measureSpec);
  NsubjettinessRatio nSub32(3, 2, axisMode, measureSpec);
  NsubjettinessRatio nSub43(4, 3, axisMode, measureSpec);
  NsubjettinessRatio nSub54(5, 4, axisMode, measureSpec);

  return {
      nSub1(pjet),  nSub2(pjet),  nSub3(pjet),  nSub4(pjet),  nSub5(pjet),
      nSub21(pjet), nSub32(pjet), nSub43(pjet), nSub54(pjet),
  };
}

void find_grandma(TClonesArray *branchGenParticles, set<int> &grandmas) {
  grandmas.clear();
  for (int i = 0; i < branchGenParticles->GetEntriesFast(); i++) {
    GenParticle *parti = (GenParticle *)branchGenParticles->At(i);
    assert(parti != NULL);
    if (parti->M1 == -1 and parti->M2 == -1) {
      // for each particles without any mother, find all the granddaughters-
      // that is the true jet
      grandmas.insert(i);
    }
  }
}

#define PARENT_UNASSIGNED -5

// TODO : try starting at the leaf nodes and going up till the grandma nodes
// (take M1 parent every time) this has a problem if different grandmas can
// create one leaf node- check if this even happens tho but how should we even
// deal with the case that a particle could have been from two different
// grandmas? ideally we would want one particle to exist in only one cluster.
void find_true_clusters(TClonesArray *branchGenParticles,
                        vector<PseudoJet> &clusters, set<int> &grandmas) {
  clusters.clear();
  clusters.reserve(grandmas.size());
  vector<set<int>> all_clusts(grandmas.size());

  // cout << "genparticle num: " << branchGenParticles->GetEntriesFast() <<
  // endl;
  vector<int> childof(branchGenParticles->GetEntriesFast());
  fill(childof.begin(), childof.end(), PARENT_UNASSIGNED);

  stack<int> doable;

  int grandma_idx = 0;
  for (const int &grandma : grandmas) {
    childof[grandma] = grandma_idx++;
    doable.push(grandma);
  }
  long long ctr = 0;
  while (not doable.empty()) {
    ctr++;

    int current;
    current = doable.top();
    doable.pop();
    assert(childof[current] != PARENT_UNASSIGNED);

    // cout << "::" << current << "\t";
    GenParticle *parti = (GenParticle *)branchGenParticles->At(current);
    assert(parti != NULL);

    vector<int> children;
    if (parti->D1 != -1)
      children.push_back(parti->D1);
    if (parti->D2 != -1 && parti->D2 != parti->D1)
      children.push_back(parti->D2);

    for (int &child : children) {
      if (childof[child] != childof[current]) {
        //  if childof[child] is PARENT_UNASSIGNED, then we are seeing it for
        //  the first time. this can occur in case :
        //    a)  it has only one parent.
        //    b)  it has 2 parents, current is the first parent. the other has
        //        not been seen yet.
        //  either way, we hit the condition as childof[current] is asserted to
        //  be assigned.
        //
        //  if childof[child] is assigned, then we hit the condition only when
        //  childof[child] is not equal to childof[current]; then we are assured
        //  that all the grandchildren have already been assigned the first
        //  parent's value. so, now, we assign them the second parent's value.
        doable.push(child);
        childof[child] = childof[current];
      }
    }

    assert(childof[current] != PARENT_UNASSIGNED);
    if (children.empty()) {
      all_clusts[childof[current]].insert(current);
    }
  }
  // NB : error in the root file: for some particles, ->D1->M1 and ->D1->M2
  // neither of them are the particle itself.
  // Dont know why this happens.
  // TODO : see why this happens

  vector<vector<PseudoJet>> all_clusts_pj(grandmas.size());
  for (int a = 0; a < all_clusts.size(); a++) {
    for (const int &element : all_clusts[a]) {
      GenParticle *parti = (GenParticle *)branchGenParticles->At(element);
      assert(parti != NULL);
      PseudoJet pj_content =
          PseudoJet(parti->Px, parti->Py, parti->Pz, parti->E);
      all_clusts_pj[a].push_back(pj_content);
    }
  }

  for (vector<PseudoJet> &clust : all_clusts_pj) {
    clusters.push_back(join(clust));
  }
}

int main(int argc, char *argv[]) {
  ofstream fout;
  run_preamble(argc, argv, fout);

  for (auto &event_title : event_types) {
    TChain chain("Delphes");
    chain.Add((rootdir + event_title + file_suffix).c_str());
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

    cout << "Processing " << event_title << endl;
    cout << "There are " << treeReader->GetEntries() << " events." << endl;

    TClonesArray *genparticle = treeReader->UseBranch("Particle");

    double ptjmin = 0; // no min pt- put cuts during analysis

    for (Long64_t entry = 0; entry < treeReader->GetEntries(); entry++) {
      treeReader->ReadEntry(entry);
      set<int> grandmas;
      find_grandma(genparticle, grandmas);
      vector<PseudoJet> true_clusters;
      find_true_clusters(genparticle, true_clusters, grandmas);

      for (PseudoJet &pjet : true_clusters) {

        double pt = pjet.pt();
        if (pt < ptjmin) {
          // cout << pt << endl;
          continue;
        }

        double mass = pjet.m();
        // double charge = pjet->Charge; // TODO : add charge ASAP
        double energy = pjet.e();
        double px = pjet.px();
        double py = pjet.py();
        double pz = pjet.pz();

        double girth = calculateGirth(pjet);

        vector<double> highero = calculateHigherMoments(pjet, girth);
        double variance = highero[0];
        double skew = highero[1];
        double kurt = highero[2];
        highero = calculateHigherMoments(pjet, girth, true);
        double var_fake = highero[0];
        double skew_fake = highero[1];
        double kurt_fake = highero[2];

        vector<double> taus = calculateNsub(pjet, 1.0);
        double tau1 = taus[0];
        double tau2 = taus[1];
        double tau3 = taus[2];
        double tau4 = taus[3];
        double tau5 = taus[4];

        double tau21 = isnan(taus[5]) ? 0 : taus[5];
        double tau32 = isnan(taus[6]) ? 0 : taus[6];
        double tau43 = isnan(taus[7]) ? 0 : taus[7];
        double tau54 = isnan(taus[8]) ? 0 : taus[8];

        double angularity_2 = calculateAngularity(pjet, 2);
        double angularity_3 = calculateAngularity(pjet, 3);
        double angularity_4 = calculateAngularity(pjet, 4);
        double angularity_5 = calculateAngularity(pjet, 5);
        double angularity_9 = calculateAngularity(pjet, 9);

        // TODO : add jet profile for some radius, planar flow, jet pull,
        // dipolarity (especially color flow)
        // TODO : add jet trajectory displacement information as per the parT
        // paper

        fout << event_title
             << ","
             // << parent_pid << ","
             << mass
             << ","
             // << charge << ","
             << pt << "," << energy << "," << px << "," << py << "," << pz
             << "," << girth << "," << variance << "," << skew << "," << kurt
             << "," << var_fake << "," << skew_fake << "," << kurt_fake << ","
             << tau1 << "," << tau2 << "," << tau3 << "," << tau4 << "," << tau5
             << "," << tau21 << "," << tau32 << "," << tau43 << "," << tau54
             << "," << angularity_2 << "," << angularity_3 << ","
             << angularity_4 << "," << angularity_5 << "," << angularity_9
             << "\n"; // newline but dont flush
      }
      if (entry % 1000 == 0) {
        fout.flush();
        cout << ".";
        cout.flush();
      }
    }
    cout << endl;
  }

  return 0;
}
