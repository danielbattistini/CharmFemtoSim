#if !defined(__CINT__) || defined(__MAKECINT__)

#include <array>
#include <string>
#include <vector>
#include <map>
#include <deque>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TSpline.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"
#include "ALICE3/Core/DelphesO2TrackSmearer.h"

#endif

using namespace Pythia8;

namespace
{
    enum tunes
    {
        kMonash = 0,
        kCRMode0,
        kCRMode2,
        kCRMode3
    };

    enum processes
    {
        kSoftQCD = 0,
        kHardQCD,
        kNonDiffractive
    };

    std::array<int, 2> DmesonPDG{411, 421}; // D+, D0
    std::array<int, 2> DstarPDG{413, 423}; // D*+, D*0
}

//__________________________________________________________________________________________________
void ComputeKstarSmearing(int nEvents=100000, int tune=kCRMode2, int process=kSoftQCD, float energy=13000, int B = 10,  int seed=42, std::string outFileNameRoot="AnalysisResults.root");
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2);

//__________________________________________________________________________________________________
void ComputeKstarSmearing(int nEvents, int tune, int process, float energy, int B, int seed, std::string outFileNameRoot)
{
    gRandom->SetSeed(seed);
    //__________________________________________________________
    // create and configure pythia generator

    Pythia pythia;
    if(process == kSoftQCD)
    {
        pythia.readString("SoftQCD:all = on");
    }
    else if(process == kNonDiffractive)
    {
        pythia.readString("SoftQCD:nonDiffractive = on");
    }
    else if(process == kHardQCD)
    {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // set tune
    if(tune == kMonash)
    {
        pythia.readString(Form("Tune:pp = 14"));
    }
    else if(tune == kCRMode0)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode2)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode3)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }

    // keep only interesting decays, to be reweighted a posteriori
    pythia.readString("421:onMode = off");
    pythia.readString("411:onMode = off");
    pythia.readString("413:onMode = off");
    pythia.readString("423:onMode = off");
    pythia.readString("421:onIfMatch = 211 321");
    pythia.readString("411:onIfMatch = 211 211 321");
    pythia.readString("413:onIfMatch = 211 421");
    pythia.readString("423:onIfMatch = 22 421"); // for simplicity, let's consider only D0* --> D0 gamma


    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    // Load Look-up tables for momentum smearing
    std::map<int, o2::delphes::TrackSmearer *> luts = {
        {211, new o2::delphes::TrackSmearer()},
        {321, new o2::delphes::TrackSmearer()},
    };

    // export CHARMFEMTOSIM=/path/to/CharmFemtoSim 
    luts[211]->loadTable(211, Form("%s/lut/lutCovm.pi.%dkG.rmin20.geometry_v2.dat", std::getenv("CHARMFEMTOSIM"), B));
    luts[321]->loadTable(321, Form("%s/lut/lutCovm.ka.%dkG.rmin20.geometry_v2.dat", std::getenv("CHARMFEMTOSIM"), B));

    //__________________________________________________________
    // define outputs
    std::map<int, std::map<int, std::map<std::string, TH2F*>>> hResoSE; // k* resolution
    for(auto &pdgDstar: DstarPDG)
    {
        for(auto &pdgDmeson: DmesonPDG)
        {
            hResoSE[pdgDstar][pdgDmeson]["part"] = new TH2F(Form("hResoSE_%d_%d", pdgDstar, pdgDmeson), "#it{k}* resolution;#it{k}*_{true} (GeV/#it{c});#it{k}*_{reco} (GeV/#it{c})", 250, 0, 0.5, 250, 0, 0.5);
            hResoSE[pdgDstar][pdgDmeson]["antipart"] = new TH2F(Form("hResoSE_%d_%d", pdgDstar, -pdgDmeson), "#it{k}* resolution;#it{k}*_{true} (GeV/#it{c});#it{k}*_{reco} (GeV/#it{c})", 250, 0, 0.5, 250, 0, 0.5);
        }
    }

    //__________________________________________________________
    // perform the simulation
    std::vector<ROOT::Math::PxPyPzMVector> partDmeson{};
    std::vector<ROOT::Math::PxPyPzMVector> partDstar{};
    std::vector<ROOT::Math::PxPyPzMVector> smearedPartDmeson{};
    std::vector<ROOT::Math::PxPyPzMVector> smearedPartDstar{};
    std::vector<int> pdgDmeson{};
    std::vector<int> pdgDstar{};
    std::vector<int> motherDmeson{};
    std::vector<int> idxDstar{};
    std::vector<float> yDmeson{};
    std::vector<float> yDstar{};

    auto fDecay = new TF1("fDecay", "exp(-x/0.1)", 0., 1000.);

    for (int iEvent=0; iEvent<nEvents; iEvent++)
    {
        int pdgD = 421;
        double massD = TDatabasePDG::Instance()->GetParticle(pdgD)->Mass();
        double ptD = 1 + gRandom->Rndm();
        double yD = 4 * (2 * gRandom->Rndm() - 1);
        double phiD = gRandom->Rndm() * 2 * TMath::Pi();
        double pxD = ptD * TMath::Cos(phiD);
        double pyD = ptD * TMath::Sin(phiD);
        double mt = TMath::Sqrt(massD * massD + ptD * ptD);
        double pzD = mt * TMath::SinH(yD);
        double pD = TMath::Sqrt(ptD * ptD + pzD * pzD);
        double ED = TMath::Sqrt(massD * massD + pD * pD);

        int DstarPdg = 413;
        double massDstar = TDatabasePDG::Instance()->GetParticle(DstarPdg)->Mass();
        double ptDstar = 1 + gRandom->Rndm();
        double y_Dstar = 4 * (2 * gRandom->Rndm() - 1);
        double phiDstar = gRandom->Rndm() * 2 * TMath::Pi();
        double pxDstar = ptDstar * TMath::Cos(phiDstar);
        double pyDstar = ptDstar * TMath::Sin(phiDstar);
        double mtDstar = TMath::Sqrt(massDstar * massDstar + ptDstar * ptDstar);
        double pzDstar = mtDstar * TMath::SinH(y_Dstar);
        double pDstar = TMath::Sqrt(ptDstar * ptDstar + pzDstar * pzDstar);
        double EDstar = TMath::Sqrt(massDstar * massDstar + pDstar * pDstar);

        // D
        Particle D;
        D.id(pdgD);
        D.status(81);
        D.m(massD);
        D.xProd(0.);
        D.yProd(0.);
        D.zProd(0.);
        D.tProd(0.);
        D.e(ED);
        D.px(pxD);
        D.py(pyD);
        D.pz(pzD);
        D.tau(fDecay->GetRandom());

        // Dstar
        Particle Dstar;
        Dstar.id(DstarPdg);
        Dstar.status(81);
        Dstar.m(massDstar);
        Dstar.xProd(0.);
        Dstar.yProd(0.);
        Dstar.zProd(0.);
        Dstar.tProd(0.);
        Dstar.e(EDstar);
        Dstar.px(pxDstar);
        Dstar.py(pyDstar);
        Dstar.pz(pzDstar);
        Dstar.tau(fDecay->GetRandom());

        pythia.event.reset();
        pythia.event.append(D);
        pythia.event.append(Dstar);
        int idPartD = pythia.event[1].id();
        int idPartDstar = pythia.event[2].id();
        pythia.particleData.mayDecay(idPartD, true);
        pythia.particleData.mayDecay(idPartDstar, true);
        pythia.moreDecays();
        
        for(int iPart=1; iPart<pythia.event.size(); iPart++)
        {
            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);
            bool isDmeson = std::find(DmesonPDG.begin(), DmesonPDG.end(), absPdg) != DmesonPDG.end();
            bool isDstar = std::find(DstarPDG.begin(), DstarPDG.end(), absPdg) != DstarPDG.end();

            // if (!isDstar && !isDmeson)
            //     continue;
            // if (isDmeson && pythia.event[iPart].pT() < 1)
            //     continue;

            auto dauList = pythia.event[iPart].daughterList();
            std::vector<double> ptDau{}, etaDau{}, pdgDau{}, eDau{};

            if(isDstar)
            {
                ROOT::Math::PxPyPzMVector part = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                partDstar.push_back(part);
                pdgDstar.push_back(pdg);
                idxDstar.push_back(iPart);
                yDstar.push_back(pythia.event[iPart].y());

                double nCh = 10; // todo: set nCh to proper multiplicity estimator

                // Smear the momentum of the daus
                std::vector<ROOT::Math::PxPyPzMVector> daus = {};
                for(const auto &iDau: dauList) {
                    auto absPdgDau = std::abs(pythia.event[iDau].id());

                    if(absPdgDau == 211 || absPdgDau == 321) {
                        double pt = pythia.event[iDau].pT();
                        double eta = pythia.event[iDau].eta();
                        
                        double ptRes =  luts[absPdgDau]->getAbsPtRes(absPdgDau, nCh, eta, pt);
                        double etaRes =  luts[absPdgDau]->getAbsEtaRes(absPdgDau, nCh, eta, pt);

                        double smearedPt = gRandom->Gaus(pt, ptRes);
                        double smearedEta = gRandom->Gaus(eta, etaRes);

                        // Assume that the smearing in px and py equally contribute to the one on pt
                        double smearedPx = pythia.event[iDau].px() * smearedPt / pt;
                        double smearedPy = pythia.event[iDau].py() * smearedPt / pt;
                        double smearedPz = pythia.event[iDau].pz(); // todo: smear pz

                        ROOT::Math::PxPyPzMVector dau(smearedPx, smearedPy, smearedPz, TDatabasePDG::Instance()->GetParticle(absPdgDau)->Mass());
                        daus.push_back(dau);
                    } else if (absPdgDau == 22) {
                        // todo: do some smearing of the photon energy
                        ROOT::Math::PxPyPzMVector dau(pythia.event[iDau].px(), pythia.event[iDau].py(), pythia.event[iDau].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                        daus.push_back(dau);
                    } else if (absPdgDau == 421) {
                        auto dauDzeroList = pythia.event[iDau].daughterList();
                        for(const auto &iDauDzero: dauDzeroList) {
                            int absPdgDauDzero = std::abs(pythia.event[iDauDzero].id());

                            double pt = pythia.event[iDauDzero].pT();
                            double eta = pythia.event[iDauDzero].eta();

                            double ptRes =  luts[absPdgDauDzero]->getAbsPtRes(absPdgDauDzero, nCh, eta, pt);
                            double etaRes =  luts[absPdgDauDzero]->getAbsEtaRes(absPdgDauDzero, nCh, eta, pt);

                            double smearedPt = gRandom->Gaus(pt, ptRes);
                            double smearedEta = gRandom->Gaus(eta, etaRes);

                            // Assume that the smearing in px and py equally contribute to the one on pt
                            double smearedPx = pythia.event[iDauDzero].px() * smearedPt / pt;
                            double smearedPy = pythia.event[iDauDzero].py() * smearedPt / pt;
                            double smearedPz = pythia.event[iDauDzero].pz(); // todo: smear pz

                            ROOT::Math::PxPyPzMVector dau(smearedPx, smearedPy, smearedPz, TDatabasePDG::Instance()->GetParticle(absPdgDauDzero)->Mass());
                            daus.push_back(dau);
                        }
                    }
                        else {
                        // todo: check
                        printf("Daughter %d not implemented. Exit!\n", absPdgDau);
                        exit(1);
                    }
                }

                ROOT::Math::PxPyPzMVector smearedPart(0, 0, 0, 0);

                // Reconstruct the D mesons with the smeared momentua
                for (const auto& dau: daus) {
                    smearedPart += dau;
                }

                // Reassign the mass of the D meson
                smearedPart = ROOT::Math::PxPyPzMVector(smearedPart.px(), smearedPart.py(), smearedPart.pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                smearedPartDstar.push_back(smearedPart);
            }
            else if(isDmeson)
            {
                ROOT::Math::PxPyPzMVector part = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                partDmeson.push_back(part);
                pdgDmeson.push_back(pdg);
                motherDmeson.push_back(pythia.event[iPart].mother1());
                yDmeson.push_back(pythia.event[iPart].y());

                double nCh = 10; // todo: set nCh to proper multiplicity estimator

                // Smear the momentum of the daus
                std::vector<ROOT::Math::PxPyPzMVector> daus = {};
                for(const auto &iDau: dauList) {
                    auto absPdgDau = std::abs(pythia.event[iDau].id());

                    if(absPdgDau == 211 || absPdgDau == 321) {
                        double pt = pythia.event[iDau].pT();
                        double eta = pythia.event[iDau].eta();
                        
                        double ptRes =  luts[absPdgDau]->getAbsPtRes(absPdgDau, nCh, eta, pt);
                        double etaRes =  luts[absPdgDau]->getAbsEtaRes(absPdgDau, nCh, eta, pt);

                        double smearedPt = gRandom->Gaus(pt, ptRes);
                        double smearedEta = gRandom->Gaus(eta, etaRes);

                        // Assume that the smearing in px and py equally contribute to the one on pt
                        double smearedPx = pythia.event[iDau].px() * smearedPt / pt;
                        double smearedPy = pythia.event[iDau].py() * smearedPt / pt;
                        double smearedPz = pythia.event[iDau].pz(); // todo: smear pz

                        ROOT::Math::PxPyPzMVector dau(smearedPx, smearedPy, smearedPz, TDatabasePDG::Instance()->GetParticle(absPdgDau)->Mass());
                        daus.push_back(dau);
                    } else {
                        // todo: check
                        printf("Daughter %d not implemented. Exit!\n", absPdgDau);
                        exit(1);
                    }
                }

                ROOT::Math::PxPyPzMVector smearedPart(0, 0, 0, 0);

                // Reconstruct the D mesons with the smeared momentua
                for (const auto& dau: daus) {
                    smearedPart += dau;
                }

                // Reassign the mass of the D meson
                smearedPart = ROOT::Math::PxPyPzMVector(smearedPart.px(), smearedPart.py(), smearedPart.pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                smearedPartDmeson.push_back(smearedPart);
            }
        }

        // same event
        for(size_t iDstar=0; iDstar<partDstar.size(); iDstar++)
        {
            for(size_t iDmeson=0; iDmeson<partDmeson.size(); iDmeson++)
            {
                 if(motherDmeson[iDmeson] == idxDstar[iDstar])
                    continue;
                std::string pair = pdgDstar[iDstar] * pdgDmeson[iDmeson] > 0 ? "part" : "antipart";

                double kStar = ComputeKstar(partDstar[iDstar], partDmeson[iDmeson]);
                double kStarSmeared = ComputeKstar(smearedPartDstar[iDstar], smearedPartDmeson[iDmeson]);
                // printf("%f  %f\n", kStar, kStarSmeared);
                hResoSE[std::abs(pdgDstar[iDstar])][std::abs(pdgDmeson[iDmeson])][pair]->Fill(kStar, kStarSmeared);
            }
        }

        partDstar.clear();
        smearedPartDstar.clear();
        pdgDstar.clear();
        idxDstar.clear();
        yDstar.clear();

        partDmeson.clear();
        smearedPartDmeson.clear();
        pdgDmeson.clear();
        motherDmeson.clear();
        yDmeson.clear();
    }

    // save root output file
    TFile outFile(outFileNameRoot.data(), "recreate");
    for(auto &pdgDstar: DstarPDG)
    {
        for(auto &pdgDmeson: DmesonPDG)
        {
            hResoSE[pdgDstar][pdgDmeson]["part"]->Write();
            hResoSE[pdgDstar][pdgDmeson]["antipart"]->Write();
        }
    }
    outFile.Close();
}

//__________________________________________________________________________________________________
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2)
{
    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);

    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;
    float kStar = 0.5 * trackRelK.P();
    return kStar;
}
