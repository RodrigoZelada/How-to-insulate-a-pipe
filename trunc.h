mesh3D trunc(mesh3D &Th, int region, map<int,int> &mapGlobalTetrahedra, map<int,int> &mapGlobalTriangles, map<int,int> &mapGlobalVertices, vector<int> &localToGlobalVertices){
    //I need maps global to local index
    
    int d=3;
    int nv=0;
    int nt=0;
    int ntet=0;
    double** vertices;
	int** triangles;
	int** tetrahedra;

    set<int> globalVertices;
    for (int n = 0; n<Th.ntet; n++){
        if (Th.tetrahedra[n][4] == region){
            mapGlobalTetrahedra.insert({n,ntet}); //global to local
            for (int i=0; i<4; i++){
                globalVertices.insert(Th.tetrahedra[n][i]);
            }
            ntet++;
        }
    }
    nv = globalVertices.size();
    localToGlobalVertices.reserve(nv);

    tetrahedra = new int* [ntet];
    for (int i = 0; i < ntet; i++) {
        tetrahedra[i] = new int[d + 2];
    }
    
    vertices = new double* [nv];
    for (int i = 0; i < nv; i++) {
        vertices[i] = new double[d+1];
    }

    nv=0;
    set<int>::iterator it;
    for (it = globalVertices.begin(); it != globalVertices.end(); it++){
        mapGlobalVertices.insert({*it, nv});
        localToGlobalVertices[nv] = *it;
        nv++;
    }
    ntet=0;
    int v;
    for (int n=0; n<Th.ntet; n++){
        if (Th.tetrahedra[n][4] == region){
            tetrahedra[ntet][4] = Th.tetrahedra[n][4];
            for (int i=0; i<4; i++){
                v =  mapGlobalVertices.at(Th.tetrahedra[n][i]);
                tetrahedra[ntet][i] = v;
                for (int j=0; j<4; j++){
                    vertices[v][j] = Th.vertices[Th.tetrahedra[n][i]][j];
                }
            }
            ntet++;
        }
    }

    /*for (int i=0; i<Th.nv; i++){
        for (int j=0; j<4; j++){
            if (mapGlobalVertices.count(i) > 0){
                vertices[mapGlobalVertices.at(i)][j] = Th.vertices[i][j];
            }
        }
    }*/

    for (int t=0; t<Th.nt; t++){
        if ( (mapGlobalVertices.count(Th.triangles[t][0]) > 0) && (mapGlobalVertices.count(Th.triangles[t][1]) > 0) && (mapGlobalVertices.count(Th.triangles[t][2]) > 0) ){
            mapGlobalTriangles.insert({t,nt});
            nt++;
        }
    }

    triangles = new int* [nt];
    for (int i = 0; i < nt; i++) {
        triangles[i] = new int[d + 1];
    }

    nt=0;
    for (int t=0; t<Th.nt; t++){
        if ( (mapGlobalVertices.count(Th.triangles[t][0]) > 0) && (mapGlobalVertices.count(Th.triangles[t][1]) > 0) && (mapGlobalVertices.count(Th.triangles[t][2]) > 0) ){
            triangles[nt][3] = Th.triangles[t][3];
            for (int i=0; i<3; i++){
                triangles[nt][i] = mapGlobalVertices.at(Th.triangles[t][i]);
            }
            nt++;
        }
    }

    mesh3D Th1(d,nv,nt,ntet,vertices,triangles,tetrahedra);
    return Th1;
}

vector<double> truncP1(vector<double> kappa, map<int,int> &mapGlobalVertices){
    int N=mapGlobalVertices.size();
    vector<double> uh(N);
    int n = kappa.size();
    for (int i=0;i<n;i++){
        if (mapGlobalVertices.count(i) > 0){
            uh[mapGlobalVertices.at(i)] = kappa[i];
        }
    }
    return uh;
}