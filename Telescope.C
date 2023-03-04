void Telescope(bool vis = true){

    TGeoManager *geom = new TGeoManager("telescope", "Telescope geometry");

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);

    //   //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
 
    //--- define the transformations
    TGeoVolume *top = geom->MakeTube("TOP", Vacuum, 0, 50., 100.); // rmin, rmax, mid height
    geom->SetTopVolume(top);

    TGeoTranslation *tr1 = new TGeoTranslation(0., 0., 50.5);
    TGeoTranslation *tr2 = new TGeoTranslation(0., 0., -50.5);

    TGeoVolume *scintillator1 = geom->MakeTube("scintillator1", Al, 0, 25., 0.5); // rmin, rmax, mid height
    scintillator1->SetLineColor(kBlue);
    top->AddNode(scintillator1, 1, tr1);
    TGeoVolume *scintillator2 = geom->MakeTube("scintillator2", Al, 0, 25., 0.5); // rmin, rmax, mid height
    scintillator2->SetLineColor(kBlue);
    top->AddNode(scintillator2, 1, tr2);

    geom->CloseGeometry();

    if (vis) top->Draw("ogle");
}