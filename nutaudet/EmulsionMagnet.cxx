#include "EmulsionMagnet.h"

#include "TGeoManager.h"
#include "FairRun.h"                    // for FairRun
#include "FairRuntimeDb.h"              // for FairRuntimeDb
#include "Riosfwd.h"                    // for ostream
#include "TList.h"                      // for TListIter, TList (ptr only)
#include "TObjArray.h"                  // for TObjArray
#include "TString.h"                    // for TString

#include "TGeoBBox.h"
#include "TGeoTrd1.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoTrd1.h"
#include "TGeoArb8.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoTransform.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include "ShipDetectorList.h"
#include "ShipUnit.h"
#include "ShipStack.h"

#include "TGeoUniformMagField.h"
#include "TVector3.h"
#include <stddef.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream,etc
#include <string.h>

using std::cout;
using std::endl;

using namespace ShipUnit;

EmulsionMagnet::~EmulsionMagnet()
{}

EmulsionMagnet::EmulsionMagnet():FairModule("EmulsionMagnet","")
{}

EmulsionMagnet::EmulsionMagnet(const char* name, const Double_t zC,const char* Title):FairModule(name, Title)
{
  fCenterZ = zC; 
}

void EmulsionMagnet::SetTPDesign(Int_t Design)
{
  fDesign = Design;
  Info("Chosen TP Design (0 no, 1 yes) "," %i", fDesign);
}

void EmulsionMagnet::SetGaps(Double_t Up, Double_t Down)
{
  fGapUpstream = Up;
  fGapDownstream = Down;
}

void EmulsionMagnet::SetMagnetSizes(Double_t X, Double_t Y, Double_t Z)
{
  fMagnetX=X;
  fMagnetY=Y;
  fMagnetZ=Z;
}

void EmulsionMagnet::SetMagnetColumn(Double_t ColX, Double_t ColY, Double_t ColZ)
{
  fColumnX=ColX;
  fColumnY=ColY;
  fColumnZ=ColZ;
}

void EmulsionMagnet::SetBaseDim(Double_t BaseX, Double_t BaseY, Double_t BaseZ)
{
  fBaseX = BaseX;
  fBaseY = BaseY;
  fBaseZ = BaseZ;
}


void EmulsionMagnet::SetCoilParameters(Double_t Radius, Double_t height1, Double_t height2, Double_t Distance)
{
  fCoilR = Radius;
  fCoilH1 = height1; //upper(left)
  fCoilH2 = height2; //lowe(right)
  fCoilDist = Distance;
}

void EmulsionMagnet::SetMagneticField(Double_t B)
{
  fField=B;
}

Int_t EmulsionMagnet::InitMedium(const char* name)
{
  static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
  static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
  static FairGeoMedia *media=geoFace->getMedia();
  static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();
    
  FairGeoMedium *ShipMedium=media->getMedium(name);
    
  if (!ShipMedium)
    {
      Fatal("InitMedium","Material %s not defined in media file.", name);
      return -1111;
    }
  TGeoMedium* medium=gGeoManager->GetMedium(name);
  if (medium!=NULL)
    return ShipMedium->getMediumIndex();
  return geoBuild->createMedium(ShipMedium);
}

void EmulsionMagnet::ConstructGeometry()
{
  TGeoVolume *top=gGeoManager->GetTopVolume();
    
  InitMedium("vacuum");
  TGeoMedium *vacuum =gGeoManager->GetMedium("vacuum");
    
  InitMedium("iron");
  TGeoMedium *Fe =gGeoManager->GetMedium("iron");
    
  InitMedium("CoilAluminium");
  TGeoMedium *Al  = gGeoManager->GetMedium("CoilAluminium");
    
  InitMedium("CoilCopper");
  TGeoMedium *Cu  = gGeoManager->GetMedium("CoilCopper");
    
  gGeoManager->SetVisLevel(10);



  cout<< "fDesign: "<< fDesign<<endl;
    
  if(fDesign==1)//OLD, TP
    {
      TGeoUniformMagField *magField1 = new TGeoUniformMagField(0.,-fField,0.); //magnetic field in Magnet pillars
      TGeoUniformMagField *magField2 = new TGeoUniformMagField(0.,fField,0.); //magnetic field in target

      TGeoBBox *MagnetBox = new TGeoBBox(fMagnetX/2, fMagnetY/2, fMagnetZ/2);
      TGeoVolume *MagnetVol = new TGeoVolume("Goliath",MagnetBox,vacuum);
      top->AddNode(MagnetVol,1,new TGeoTranslation(0,0,fCenterZ));
      
      //Iron basis on which the coils are placed
      TGeoBBox *Base = new TGeoBBox(fBaseX/2,fBaseY/2,fBaseZ/2);
      TGeoVolume *volBase = new TGeoVolume("volBase",Base,Fe);
      volBase->SetLineColor(kRed);
      MagnetVol->AddNode(volBase,1,new TGeoTranslation(0, fMagnetY/2 - fBaseY/2, 0)); //upper part
      MagnetVol->AddNode(volBase,2,new TGeoTranslation(0, -fMagnetY/2 + fBaseY/2, 0)); //lower part

      //Coils Description: 2 volumes must be defined being the upper coil in Cu and the lower one in Al and also heghts are different
      TGeoTube *CoilBoxU = new TGeoTube("C",0,fCoilR,fCoilH1/2);
      TGeoVolume *CoilVolUp = new TGeoVolume("CoilVolUp",CoilBoxU, Cu);
      CoilVolUp->SetLineColor(kGreen);
      TGeoTube *CoilBoxD = new TGeoTube("C",0,fCoilR,fCoilH2/2);
      TGeoVolume *CoilVolDown = new TGeoVolume("CoilVolDown",CoilBoxD, Al);
      CoilVolDown->SetLineColor(kGreen);
      
      TGeoRotation *r1 = new TGeoRotation();
      r1->SetAngles(0,90,0);
      TGeoCombiTrans tUp(0, fMagnetY/2 - fBaseY - fCoilH1/2, 0,r1);
      TGeoHMatrix *mUp = new TGeoHMatrix(tUp);
      TGeoCombiTrans tDown(0, -fMagnetY/2 + fBaseY + fCoilH2/2, 0,r1);
      TGeoHMatrix *mDown = new TGeoHMatrix(tDown);

      MagnetVol->AddNode(CoilVolUp,1,mUp);
      MagnetVol->AddNode(CoilVolDown,1,mDown);
      
      //********************* Columns ****************************
    
      //Each column is made of a longer pillar (rectangle + trapezoid) and on top a shorter pillar (rectangle + trapezoid again) 
     
    Double_t base1 = 135, base2 = 78; //basis of the trapezoid
    Double_t side1 = 33, side2 = 125, side3 = 57, side4 = 90; //Sides of the columns
    
    //***** SIDE Left Front ****
    
    //Shorter Pillar: rectangle
    TGeoBBox *LateralS1 = new TGeoBBox("LateralS1",side1/2,fCoilH1/2,base1/2);
    TGeoTranslation *tr1 = new TGeoTranslation(-fMagnetX/2 + side1/2, fMagnetY/2 - fBaseY - fCoilH1/2, -fMagnetZ/2 + base1/2);
    TGeoVolume *volLateralS1 = new TGeoVolume("volLateralS1",LateralS1,Fe);
    volLateralS1->SetLineColor(kRed);
    volLateralS1->SetField(magField1);
    MagnetVol->AddNode(volLateralS1, 1, tr1);
    
    //Shorter Pillar: trapezoid
    
    TGeoArb8 *LateralS2 = new TGeoArb8("LateralS2",fCoilH1/2);
    LateralS2->SetVertex(0, side4, 0);
    LateralS2->SetVertex(1, side1, 0);
    LateralS2->SetVertex(2, side1, base1);
    LateralS2->SetVertex(3, side4, base2);
    LateralS2->SetVertex(4, side4, 0);
    LateralS2->SetVertex(5, side1, 0);
    LateralS2->SetVertex(6, side1, base1);
    LateralS2->SetVertex(7, side4, base2);
    
    TGeoVolume *volLateralS2 = new TGeoVolume("volLateralS2",LateralS2,Fe);
    volLateralS2->SetLineColor(kRed);
    volLateralS2->SetField(magField1);
    
    TGeoRotation *r2 = new TGeoRotation();
    r2->SetAngles(0,90,0);
    TGeoCombiTrans tr3(-fMagnetX/2, fMagnetY/2 - fBaseY - fCoilH1/2, -fMagnetZ/2,r2);
    TGeoHMatrix *m3_a = new TGeoHMatrix(tr3);
    MagnetVol->AddNode(volLateralS2, 1, m3_a);

    //LOWER LATERAL SURFACE
    
    //LONGER RECTANGLE
    TGeoBBox *LateralSurface1low = new TGeoBBox("LateralSurface1low",side1/2,(fCoilDist + fCoilH2)/2,side2/2);
    TGeoVolume *volLateralSurface1low = new TGeoVolume("volLateralSurface1low",LateralSurface1low,Fe);
    volLateralSurface1low->SetLineColor(kRed);
    volLateralSurface1low->SetField(magField1);
    TGeoTranslation *tr1low = new TGeoTranslation(-fMagnetX/2 +side1/2, fMagnetY/2 - fBaseY - fCoilH1 - (fCoilDist + fCoilH2)/2, -fMagnetZ/2 + side2/2);
    MagnetVol->AddNode(volLateralSurface1low, 1, tr1low);;
    
    
    //SHORTER RECTANGLE
    TGeoBBox *LateralSurface2low = new TGeoBBox("LateralSurface2low",side3/2,(fCoilDist + fCoilH2)/2,base2/2);
    TGeoVolume *volLateralSurface2low = new TGeoVolume("volLateralSurface2low",LateralSurface2low,Fe);
    volLateralSurface2low->SetLineColor(kRed);
    TGeoTranslation *tr2low = new TGeoTranslation(-fMagnetX/2 +side1 + side3/2, fMagnetY/2 - fBaseY -fCoilH1 - (fCoilDist + fCoilH2)/2, -fMagnetZ/2 + base2/2);
    MagnetVol->AddNode(volLateralSurface2low, 1, tr2low);
    volLateralSurface2low->SetField(magField1);

    //***** SIDE Right Front ****
    
    //LONGER RECTANGLE
    TGeoTranslation *tr1_b = new TGeoTranslation(-fMagnetX/2 + side1/2, fMagnetY/2 - fBaseY - fCoilH1/2, fMagnetZ/2 - base1/2);
    TGeoVolume *volLateralS1_b = new TGeoVolume("volLateralS1_b",LateralS1,Fe);
    volLateralS1_b->SetLineColor(kRed);
    volLateralS1_b->SetField(magField1);
    MagnetVol->AddNode(volLateralS1_b, 1, tr1_b);
    
    //TRAPEZOID
    TGeoArb8 *LateralS2_b = new TGeoArb8("LateralS2_b",fCoilH1/2);
    LateralS2_b ->SetVertex(0, side4, 0);
    LateralS2_b ->SetVertex(1, side1, 0);
    LateralS2_b ->SetVertex(2, side1, base1);
    LateralS2_b ->SetVertex(3, side4, base2);
    LateralS2_b ->SetVertex(4, side4, 0);
    LateralS2_b ->SetVertex(5, side1, 0);
    LateralS2_b ->SetVertex(6, side1, base1);
    LateralS2_b ->SetVertex(7, side4, base2);
    
    TGeoVolume *volLateralS2_b = new TGeoVolume("volLateralS2_b",LateralS2_b,Fe);
    volLateralS2_b->SetLineColor(kRed);
    volLateralS2_b->SetField(magField1);
    
    TGeoRotation *r2_b = new TGeoRotation();
    r2_b->SetAngles(0,270,0);
    TGeoCombiTrans tr2_b(-fMagnetX/2 , fMagnetY/2 - fBaseY - fCoilH1/2, fMagnetZ/2,r2_b);
    TGeoHMatrix *m3_b = new TGeoHMatrix(tr2_b);
    MagnetVol->AddNode(volLateralS2_b, 1, m3_b);
    
    
    //LOWER LATERAL SURFACE
    
    //LONGER RECTANGLE
    TGeoVolume *volLateralSurface1blow = new TGeoVolume("volLateralSurface1blow",LateralSurface1low,Fe);
    volLateralSurface1blow->SetLineColor(kRed);
    volLateralSurface1blow->SetField(magField1);
    TGeoTranslation *tr1blow = new TGeoTranslation(-fMagnetX/2 +side1/2, fMagnetY/2 - fBaseY - fCoilH1 - (fCoilDist + fCoilH2)/2, fMagnetZ/2 - side2/2);
    MagnetVol->AddNode(volLateralSurface1blow, 1, tr1blow);;
    
    
    //SHORTER RECTANGLE
    TGeoVolume *volLateralSurface2blow = new TGeoVolume("volLateralSurface2blow",LateralSurface2low,Fe);
    volLateralSurface2blow->SetLineColor(kRed);
    volLateralSurface2blow->SetField(magField1);
    TGeoTranslation *tr2blow = new TGeoTranslation(-fMagnetX/2 +side1 + side3/2, fMagnetY/2 - fBaseY - fCoilH1 - (fCoilDist + fCoilH2)/2, fMagnetZ/2 - base2/2);
    MagnetVol->AddNode(volLateralSurface2blow, 1, tr2blow);
    
    
    //***** SIDE left Back ****
    
    
    //LONGER RECTANGLE
    TGeoBBox *LateralS1_d = new TGeoBBox("LateralS1_d",side1/2,(fCoilH1 + fCoilH2 + fCoilDist)/2,base1/2);
    TGeoTranslation *tr1_d = new TGeoTranslation(fMagnetX/2 - side1/2, fMagnetY/2 - fBaseY - (fCoilH1 + fCoilH2 + fCoilDist)/2, -fMagnetZ/2 + base1/2);
    TGeoVolume *volLateralS1_d = new TGeoVolume("volLateralS1_d",LateralS1_d,Fe);
    volLateralS1_d->SetLineColor(kRed);
    volLateralS1_d->SetField(magField1);
    MagnetVol->AddNode(volLateralS1_d, 1, tr1_d);
    
    //TRAPEZOID
    
    TGeoArb8 *LateralS2_d = new TGeoArb8("LateralS2_d",(fCoilH1 + fCoilH2 + fCoilDist)/2);
    LateralS2_d->SetVertex(0, side4, 0);
    LateralS2_d->SetVertex(1, side1, 0);
    LateralS2_d->SetVertex(2, side1, base1);
    LateralS2_d->SetVertex(3, side4, base2);
    LateralS2_d->SetVertex(4, side4, 0);
    LateralS2_d->SetVertex(5, side1, 0);
    LateralS2_d->SetVertex(6, side1, base1);
    LateralS2_d->SetVertex(7, side4, base2);
    
    
    TGeoVolume *volLateralS2_d = new TGeoVolume("volLateralS2_d",LateralS2_d,Fe);
    volLateralS2_d->SetLineColor(kRed);
    volLateralS2_d->SetField(magField1);
    
    TGeoRotation *r2_d = new TGeoRotation();
    r2_d->SetAngles(0,270,180);
    TGeoCombiTrans tr2_d(fMagnetX/2 , fMagnetY/2 - fBaseY - (fCoilH1 + fCoilH2 + fCoilDist)/2, -fMagnetZ/2,r2_d);
    TGeoHMatrix *m3_d = new TGeoHMatrix(tr2_d);
    MagnetVol->AddNode(volLateralS2_d, 1, m3_d);

//***** SIDE right Back ****
    
    
    //LONGER RECTANGLE
    TGeoBBox *LateralS1_c = new TGeoBBox("LateralS1_c",side1/2,(fCoilH1 + fCoilH2 + fCoilDist)/2,base1/2);
    TGeoTranslation *tr1_c = new TGeoTranslation(fMagnetX/2 - side1/2, fMagnetY/2 - fBaseY - (fCoilH1 + fCoilH2 + fCoilDist)/2, fMagnetZ/2 - base1/2);
    TGeoVolume *volLateralS1_c = new TGeoVolume("volLateralS1_c",LateralS1_c,Fe);
    volLateralS1_c->SetLineColor(kRed);
    volLateralS1_c->SetField(magField1);
    MagnetVol->AddNode(volLateralS1_c, 1, tr1_c);
    
    //TRAPEZOID
    
    TGeoArb8 *LateralS2_c = new TGeoArb8("LateralS2_c",(fCoilH1 + fCoilH2 + fCoilDist)/2);
    LateralS2_c ->SetVertex(0, side4, 0);
    LateralS2_c ->SetVertex(1, side1, 0);
    LateralS2_c ->SetVertex(2, side1, base1);
    LateralS2_c ->SetVertex(3, side4, base2);
    LateralS2_c ->SetVertex(4, side4, 0);
    LateralS2_c ->SetVertex(5, side1, 0);
    LateralS2_c ->SetVertex(6, side1, base1);
    LateralS2_c ->SetVertex(7, side4, base2);
    
    TGeoVolume *volLateralS2_c = new TGeoVolume("volLateralS2_c",LateralS2_c,Fe);
    volLateralS2_c->SetLineColor(kRed);
    volLateralS2_c->SetField(magField1);
    
    TGeoRotation *r2_c = new TGeoRotation();
    r2_c->SetAngles(0,90,180);
    TGeoCombiTrans tr2_c(fMagnetX/2 , fMagnetY/2 - fBaseY - (fCoilH1 + fCoilH2 + fCoilDist)/2, fMagnetZ/2,r2_c);
    TGeoHMatrix *m3_c = new TGeoHMatrix(tr2_c);
    MagnetVol->AddNode(volLateralS2_c, 1, m3_c);


    }

  if(fDesign==0) //NEW
    {
      TGeoUniformMagField *magField1 = new TGeoUniformMagField(-fField,0.,0.); //magnetic field in Magnet pillars
      TGeoUniformMagField *magField2 = new TGeoUniformMagField(fField,0.,0.); //magnetic field in target
      
      TGeoBBox *MagnetBox = new TGeoBBox(fMagnetX/2, fMagnetY/2, fMagnetZ/2);
      TGeoVolume *MagnetVol = new TGeoVolume("Davide",MagnetBox,vacuum);
      top->AddNode(MagnetVol,1,new TGeoTranslation(0,0,fCenterZ));
    
      //The -0.01*mm is only for drawing reasons
      TGeoBBox *LateralBox = new TGeoBBox("LB",fBaseZ/2,fBaseY/2,(fBaseX-0.01*mm)/2);
      TGeoTube *CoilBox = new TGeoTube("C",0,fCoilR,fCoilH1/2);
    
      TGeoCompositeShape *LateralSurf = new TGeoCompositeShape("LS","LB-C");
    
      TGeoVolume *CoilVol = new TGeoVolume("CoilVol",CoilBox, Cu);
      CoilVol->SetLineColor(kGreen);
      TGeoVolume *LateralSurfVol = new TGeoVolume("LateralSurfVol",LateralSurf, Fe);
      LateralSurfVol->SetLineColor(kRed);

      TGeoRotation *r1 = new TGeoRotation();
      r1->RotateY(90);
      //r1->RotateY(90);
      //r1->RotateY(90);
      r1->RegisterYourself();
      TGeoTranslation *t1r = new TGeoTranslation(-fMagnetX/2+fBaseX/2,0,0);
      TGeoTranslation *t1l = new TGeoTranslation(fMagnetX/2-fBaseX/2,0,0);

      TGeoCombiTrans *trans1r = new TGeoCombiTrans(-fMagnetX/2+fBaseX/2,0,0,r1);
      TGeoCombiTrans *trans1l = new TGeoCombiTrans(fMagnetX/2-fBaseX/2,0,0,r1);
      TGeoHMatrix *m1_r = new TGeoHMatrix("m1_r");
      *m1_r = trans1r;
      TGeoHMatrix *m1_l = new TGeoHMatrix("m1_l");
      *m1_l = trans1l;

      MagnetVol->AddNode(CoilVol,1, m1_r);
      MagnetVol->AddNode(LateralSurfVol,1,m1_r);
      MagnetVol->AddNode(CoilVol,2, m1_l);
      MagnetVol->AddNode(LateralSurfVol,2,m1_l);

      TGeoBBox *ColumnBox = new TGeoBBox(fColumnX/2, fColumnY/2, fColumnZ/2);
      TGeoVolume *ColumnVol = new TGeoVolume("ColumnVol",ColumnBox,Fe);
      ColumnVol->SetField(magField1);
      ColumnVol->SetLineColor(kRed);
      MagnetVol->AddNode(ColumnVol,1,new TGeoTranslation(0,fMagnetY/2-fColumnY/2, -fMagnetZ/2+fColumnZ/2));
      MagnetVol->AddNode(ColumnVol,2,new TGeoTranslation(0,fMagnetY/2-fColumnY/2, fMagnetZ/2-fColumnZ/2));
      MagnetVol->AddNode(ColumnVol,3,new TGeoTranslation(0,-fMagnetY/2+fColumnY/2, -fMagnetZ/2+fColumnZ/2));
      MagnetVol->AddNode(ColumnVol,4,new TGeoTranslation(0,-fMagnetY/2+fColumnY/2, fMagnetZ/2-fColumnZ/2));
    
      TGeoBBox *BaseBox = new TGeoBBox(fCoilDist/2,fColumnY/2, fBaseZ/2);
      TGeoVolume *BaseVol = new TGeoVolume("BaseVol",BaseBox,Fe);
      BaseVol->SetLineColor(kRed);
      MagnetVol->AddNode(BaseVol,1, new TGeoTranslation(0,-fMagnetY/2+fColumnY/2,0));
    }
}

ClassImp(EmulsionMagnet)