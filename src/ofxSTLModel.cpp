// adapted from the STL class by Marius Watz (part of unlekkerlib from processing)
// taken from http://workshop.evolutionzone.com/unlekkerlib/ on 06/12/09
// converted to C++ by Marek Bereza
// slightly tweaked by George Profenza

#include "ofxSTLModel.h"

	
	
ofxSTLModel::ofxSTLModel() {
	mult = 1;
}

bool ofxSTLModel::validSTL(string _path) {
    
    // Load the file
    ofFile inFile;
    if (!inFile.open(_path, ofFile::ReadOnly)) {
        inFile.close();
        return false;
    }
    
    // Check for a valid ASCII file
    char* bytes = new char[6];
    inFile.read(bytes, 6);
    std::string beg( reinterpret_cast< char const* >(bytes) ) ;
    if (beg == "solid ") {
        inFile.close();
        return true;
    }
    
    // Check for a valid binary file
    inFile.seekg(80);
    unsigned int nFacets;
    inFile.read((char*)(&nFacets), 4);
    bool bBinaryOK = (84 + nFacets * 50) == inFile.getSize();
    inFile.close();
    return bBinaryOK;
}


/**
 * Draws the object.
 */
void ofxSTLModel::draw() {
	vboMesh.draw();
}

/**
 * Calculates the bounding box of the object. 
 */
	
void ofxSTLModel::calcBounds() {
	minx=10000;
	maxx=-10000;
	miny=10000;
	maxy=-10000;
	minz=10000;
	maxz=-10000;
	
	int id=0;
	for(int i=0; i<triangles.size(); i++) {
		id=3;
		for(int j=0; j<3; j++) {
			if(triangles[i].v[id]<minx) minx=triangles[i].v[id];
			else if(triangles[i].v[id]>maxx) maxx=triangles[i].v[id];
			id++;
			
			if(triangles[i].v[id]<miny) miny=triangles[i].v[id];
			else if(triangles[i].v[id]>maxy) maxy=triangles[i].v[id];
			id++;
			
			if(triangles[i].v[id]<minz) minz=triangles[i].v[id];
			else if(triangles[i].v[id]>maxz) maxz=triangles[i].v[id];
			id++;
			
		}
	}
}
	
/**
 * Centers the object around the world origin. 
 */

void ofxSTLModel::center() {
	calcBounds();
	float tx=(minx+maxx)/2;
	float ty=(miny+maxy)/2;
	float tz=(minz+maxz)/2;
	for(int i=0; i<triangles.size(); i++) triangles[i].translate(-tx,-ty,-tz);
}

	
/**
 * Normalizes the object to a absolute scale 
 */

void ofxSTLModel::normalize(float m) {
	calcBounds();
	
	float max=maxx-minx;
	if(maxy-miny>max) max=maxy-miny;
	if(maxz-minz>max) max=maxz-minz;
	
	float tx=m/max;
	float ty=m/max;
	float tz=m/max;
	for(int i=0; i<triangles.size(); i++) triangles[i].scale(tx,ty,tz);
}
	
	
	/////////////////////////////////////////////
	// FUNCTIONS FOR STL INPUT
	
void ofxSTLModel::read(string path) {
	
	char header[80];
	char buf[50]; // 12 32bit floats + 2 bytes = 50 bytes per triangle
	int num = 0;
	
	ifstream file(ofToDataPath(path).c_str(), ios::in|ios::binary);
	
	printf("Reading %s\n", path.c_str());
	
	file.read(header, 80);
	file.read((char*)&num, 4);
	
	printf("%d facets\n", num);
	triangles.clear();
	
	for(int i=0; i<num; i++) {
		file.read(buf, 50);
		triangles.push_back(STLFace());
		triangles[i].parseFace(buf);
		//if(i%1000==0) printf("%d triangles read.", i);
	}
	file.close();
    
    vboMesh.setMode(OF_PRIMITIVE_TRIANGLES);
	
	for(int i=0; i< triangles.size(); i++) {
		vboMesh.addNormal(ofVec3f(triangles[i].v[0], triangles[i].v[1], triangles[i].v[2]));
//        vboMesh.addNormal(ofVec3f(triangles[i].v[0], triangles[i].v[1], triangles[i].v[2]));
//        vboMesh.addNormal(ofVec3f(triangles[i].v[0], triangles[i].v[1], triangles[i].v[2]));
		vboMesh.addVertex(ofVec3f(triangles[i].v[3], triangles[i].v[4], triangles[i].v[5]));
		vboMesh.addVertex(ofVec3f(triangles[i].v[6], triangles[i].v[7], triangles[i].v[8]));
		vboMesh.addVertex(ofVec3f(triangles[i].v[9], triangles[i].v[10], triangles[i].v[11]));
//        vboMesh.addTriangle(vboMesh.getVertices().size()-3, vboMesh.getVertices().size()-2, vboMesh.getVertices().size()-1);
	}
}


void ofxSTLModel::readToMesh(string path, ofVboMesh* mesh, int nEveryTriangles) {
    
    bool bBinary = true;
    ofFile inFile;
    inFile.open(path, ofFile::ReadOnly);
    char* bytes = new char[6];
    inFile.read(bytes, 6);
    std::string beg( reinterpret_cast< char const* >(bytes) ) ;
    if (beg == "solid ") {
        inFile.close();
        bBinary = false;
    }
    
    if (bBinary) {
        cout << "Reading Binary File" << endl;
        readToMeshBinary(path, mesh, nEveryTriangles);
    } else {
        cout << "Reading ASCII File" << endl;
        readToMeshAscii(path, mesh, nEveryTriangles);
    }
}

void ofxSTLModel::readToMeshAscii(string path, ofVboMesh* mesh, int nEveryTriangles) {
    
    // Clear the provided mesh
    mesh->clear();
    mesh->setMode(OF_PRIMITIVE_TRIANGLES);
    mesh->disableColors();
    mesh->disableTextures();
    
    
    ifstream file(ofToDataPath(path).c_str(), ios::in);
    std::string line;
    string vertID = "vertex ";
    vertices.clear();
    while(std::getline(file, line)) {
        
        int loc = line.find(vertID);
        if (loc != -1) {
            
            // Add a new vertex
            vertices.push_back(ofVec3f(0, 0, 0));
            
            // Find the start of the numbers
            int pre = loc + vertID.length();
            int len = line.length() - pre;
            line = line.substr(pre);
            
            // Get all of the coordinates
            loc = line.find(' ');
            vertices.back().x = ofToFloat(line.substr(0, loc));
            line = line.substr(loc+1);
            loc = line.find(' ');
            vertices.back().y = ofToFloat(line.substr(0, loc));
            vertices.back().z = ofToFloat(line.substr(loc));
        }
    }
    
    if (vertices.size() % 3 != 0) {
        cout << "Wrong number of vertices in ASCII STL file." << endl;
    }
    int nFacets = vertices.size()/3;
    int nVerts = nFacets*3;
    int nNorms = nFacets;
    float* verts = (float*)(&(vertices[0]));
    
    // Normalize all facets to a spherical bounding box
    float radius;
    ofVec3f center;
    bounding_sphere(radius, (float*)&(center[0]), verts, nFacets*3);
    
    ofVec3f m, M;
    m = center - ofVec3f(radius, radius, radius);
    M = center + ofVec3f(radius, radius, radius);
    
    for (int i = 0; i < nFacets*3; i++) {
        
        ofVec3f* pt = &(vertices[i]);
        pt->x = ofMap(pt->x, m.x, M.x, minVal, maxVal, true);
        pt->y = ofMap(pt->y, m.y, M.y, minVal, maxVal, true);
        pt->z = ofMap(pt->z, m.z, M.z, minVal, maxVal, true);
    }
    
    // Save every nth triangle to a mesh for viewing
    for (int i = 0; i < nFacets; i += nEveryTriangles) {
        
        mesh->addVertex(vertices[i*3+0]);
        mesh->addVertex(vertices[i*3+1]);
        mesh->addVertex(vertices[i*3+2]);
    }
    
    // Close the file
    file.close();
}



void ofxSTLModel::readToMeshBinary(string path, ofVboMesh* mesh, int nEveryTriangles) {
    
    // Clear the provided mesh
    mesh->clear();
    mesh->setMode(OF_PRIMITIVE_TRIANGLES);
//    mesh->disableColors();
//    mesh->disableTextures();
    
    // Data storage
    char header[80];
    char facetData[50]; // 12 32bit floats + 2 bytes = 50 bytes per triangle
    int nFacets = 0;
    float data[12];
    
    // Open File and begin reading
    ifstream file(ofToDataPath(path).c_str(), ios::in|ios::binary);
    file.read(header, 80);
    file.read((char*)&nFacets, 4);
    
    int nVerts = nFacets*3;
    int nNorms = nFacets;
    
    // Iterate through all facets
    vertices.clear();
    vertices.resize(nVerts);
    float* verts = (float*)(&(vertices[0]));
    normals.clear();
    normals.resize(nFacets);
    float* norms = (float*)(&(normals[0]));
    for (int i = 0; i < nFacets; i++) {
        
        // Read the data of this facet
        file.read(facetData, 50);
        memcpy(data, facetData, 12*sizeof(float));
        
        // Save the vert data
        memcpy((char*)verts+i*9*sizeof(float), (char*)data+3*sizeof(float), 9*sizeof(float));
        
        // Save the normal data
        memcpy((char*)norms+i*3*sizeof(float), (char*)data, 3*sizeof(float));
    }
    
    // Normalize all facets to a spherical bounding box
    float radius;
    ofVec3f center;
    bounding_sphere(radius, (float*)&(center[0]), verts, nFacets*3);
    
    ofVec3f m, M;
    m = center - ofVec3f(radius, radius, radius);
    M = center + ofVec3f(radius, radius, radius);
    
    for (int i = 0; i < nFacets*3; i++) {
        
        ofVec3f* pt = &(vertices[i]);
        pt->x = ofMap(pt->x, m.x, M.x, minVal, maxVal, true);
        pt->y = ofMap(pt->y, m.y, M.y, minVal, maxVal, true);
        pt->z = ofMap(pt->z, m.z, M.z, minVal, maxVal, true);
    }
    
    // Save every nth triangle to a mesh for viewing
  
    for (int i = 0; i < nFacets; i += nEveryTriangles) {
    
        // This loads correctly
//        cout << "Face " << i << endl;
//        cout << "\t" << "vert\t" << vertices[i*3+0] << endl;
//        cout << "\t" << "vert\t" << vertices[i*3+1] << endl;
//        cout << "\t" << "vert\t" << vertices[i*3+2] << endl;
//        cout << "\t" << "norm\t" << normals[i] << endl;
        
        
        
        mesh->addVertex(vertices[i*3+0]);
        mesh->addVertex(vertices[i*3+1]);
        mesh->addVertex(vertices[i*3+2]);
        
//        mesh->addNormal(getSTLNormal(vertices[i*3+0],
//                                     vertices[i*3+1],
//                                     vertices[i*3+2]));
//        mesh->addNormal(getSTLNormal(vertices[i*3+0],
//                                     vertices[i*3+1],
//                                     vertices[i*3+2]));
//        mesh->addNormal(getSTLNormal(vertices[i*3+0],
//                                     vertices[i*3+1],
//                                     vertices[i*3+2]));
        mesh->addNormal(normals[i]);
        mesh->addNormal(normals[i]);
        mesh->addNormal(normals[i]);
        
        mesh->addColor(ofColor(255));
        mesh->addColor(ofColor(255));
        mesh->addColor(ofColor(255));
    }
//    mesh->addVertex(ofVec3f(0, 0, 0));
//    mesh->addVertex(ofVec3f(100, 0, 0));
//    mesh->addVertex(ofVec3f(100, 100, 0));
//    mesh->addNormal(ofVec3f(0, 0, -1));
//    mesh->addNormal(ofVec3f(0, 0, -1));
//    mesh->addNormal(ofVec3f(0, 0, -1));
//    mesh->addColor(ofColor(255));
//    mesh->addColor(ofColor(255));
//    mesh->addColor(ofColor(255));
    
    // Close the file
    file.close();
}
	
	
	
/////////////////////////////////////////////
// FUNCTIONS FOR RAW STL OUTPUT


//--------------------------------------------------------------
void ofxSTLModel::bounding_sphere(float& radius, float* center, const float* V, const int n)
{
    int d = 3; // 3D mini-ball
    radius = center[0] = center[1] = center[2] = 0;
    if (n < 2) return;
    
    // mini-ball
    const float** ap = new const float*[n];
    for (int i = 0; i < n; ++i) { ap[i] = V + d * i; }
    typedef const float** PointIterator;
    typedef const float* CoordIterator;
    Miniball::Miniball <
    Miniball::CoordAccessor < PointIterator, CoordIterator >>
    miniball(d, ap, ap + n);
    
    // get result
    if (miniball.is_valid())
    {
        const float* cnt = miniball.center();
        for (int i = 0; i < d; ++i) {
            center[i] = cnt[i];
        }
        radius = sqrtf(miniball.squared_radius() + 1.0e-20f);
    }
    else
    {
        // the miniball might failed sometimes
        // if so, just calculate the bounding box
        
        float bbmin[3] = { V[0], V[1], V[2] };
        float bbmax[3] = { V[0], V[1], V[2] };
        for (int i = 1; i < n; ++i)
        {
            int i3 = i * 3;
            for (int j = 0; j < d; ++j)
            {
                float tmp = V[i3 + j];
                if (tmp < bbmin[j]) bbmin[j] = tmp;
                if (tmp > bbmax[j]) bbmax[j] = tmp;
            }
        }
        
        float width[3];
        for (int j = 0; j < d; ++j)
        {
            width[j] = (bbmax[j] - bbmin[j]) / 2.0f;
            center[j] = (bbmax[j] + bbmin[j]) / 2.0f;
        }
        
        radius = width[0];
        if (width[1] > radius) radius = width[1];
        if (width[2] > radius) radius = width[2];
    }
    
    // release
    delete[] ap;
}


void ofxSTLModel::write(string path) {
	int num = triangles.size();
	ofstream file(path.c_str(), ofstream::binary);
	
	char header[80];
	
	memset(header, 0, 80);
	file.write(header, 80);
	file.write((char*)&num, 4);
	char buf[50];
	memset(header, 0, 50);
	for(int i = 0; i < triangles.size(); i++) {
		memcpy(buf, triangles[i].getData(), 48);
		file.write(buf, 50);
	}
	
	file.close();
}


void ofxSTLModel::writeFromMesh(string path, float angle, ofVec3f pivot, ofVec3f axis) {
    
    int nVerts = vertices.size();
    unsigned int nFacets = nVerts/3;
    int nNorms = nFacets;
    
    // Create an output file
    ofstream file(path.c_str(), ofstream::binary);
    
    // Write the header to file
    char header[80];
    memset(header, 0, 80);
    file.write(header, 80);
    file.write((char*)&nFacets, 4);
    
    // Iterate over all vertices, transforming them
    char buf[50];
    memset(buf, 0, 50);
    for (int i = 0; i < nFacets; i++) {
        
        // Rotate each vertex
        for (int j = 0; j < 3; j++) {
            vertices[i*3+j] = vertices[i*3+j].getRotated(angle, pivot, axis);
        }
        
        // Calculate the normal
        ofVec3f normal = getSTLNormal((float*)(&(vertices[i*3])));
        
        // Save these values to file
        memcpy(buf, (char*)(&(normal[0])), 3*sizeof(float));
        memcpy(buf+3*sizeof(float), (char*)(&(vertices[i*3])), 9*sizeof(float));
        file.write(buf, 50);
    }
    file.close();
}

void ofxSTLModel::writeFromMesh(string path, float scale, float ax, float ay, float az, ofVec3f translation) {
    
    int nVerts = vertices.size();
    unsigned int nFacets = nVerts/3;
    int nNorms = nFacets;
    
    // Create an output file
    ofstream file(path.c_str(), ofstream::binary);
    
    // Write the header to file
    char header[80];
    memset(header, 0, 80);
    file.write(header, 80);
    file.write((char*)&nFacets, 4);
    
    // Iterate over all vertices, transforming them
    char buf[50];
    memset(buf, 0, 50);
    for (int i = 0; i < nFacets; i++) {
        
        // Rotate each vertex
        for (int j = 0; j < 3; j++) {
            vertices[i*3+j] = (vertices[i*3+j]*scale).getRotated(ax, ay, az)+translation;
        }
        
        // Calculate the normal
        ofVec3f normal = getSTLNormal((float*)(&(vertices[i*3])));
        
        // Save these values to file
        memcpy(buf, (char*)(&(normal[0])), 3*sizeof(float));
        memcpy(buf+3*sizeof(float), (char*)(&(vertices[i*3])), 9*sizeof(float));
        file.write(buf, 50);
    }
    file.close();
}

ofVec3f ofxSTLModel::getSTLNormal(float* v) {
    
    float x,y,z,vx1,vy1,vz1,vx2,vy2,vz2;
    vx1=v[3]-v[0];
    vy1=v[4]-v[1];
    vz1=v[5]-v[2];
    vx2=v[6]-v[0];
    vy2=v[7]-v[1];
    vz2=v[8]-v[2];
    
    x=vy1*vz2-vy2*vz1;
    y=vz1*vx2-vz2*vx1;
    z=vx1*vy2-vx2*vy1;
    
    float l=(float)sqrt(x*x+y*y+z*z);
    x/=l;
    y/=l;
    z/=l;
    return ofVec3f(x, y, z);
}

ofVec3f ofxSTLModel::getSTLNormal(ofVec3f v1, ofVec3f v2, ofVec3f v3) {
    
    return (v2-v1).getCrossed(v3-v1).normalize();
    
    
    
    float* tmp = new float[9];
    tmp[0] = v1.x;
    tmp[1] = v1.y;
    tmp[2] = v1.z;
    tmp[3] = v2.x;
    tmp[4] = v2.y;
    tmp[5] = v2.z;
    tmp[6] = v3.x;
    tmp[7] = v3.y;
    tmp[8] = v3.z;
    
    return getSTLNormal(tmp);
}


//float x,y,z,vx1,vy1,vz1,vx2,vy2,vz2;
//vx1=v[6]-v[3];
//vy1=v[7]-v[4];
//vz1=v[8]-v[5];
//vx2=v[9]-v[3];
//vy2=v[10]-v[4];
//vz2=v[11]-v[5];
//
//x=vy1*vz2-vy2*vz1;
//y=vz1*vx2-vz2*vx1;
//z=vx1*vy2-vx2*vy1;
//
//float l=(float)sqrt(x*x+y*y+z*z);
//x/=l;
//y/=l;
//z/=l;
//v[0]=x;
//v[1]=y;
//v[2]=z;

		
	
void ofxSTLModel::addTriangle(float nX, float nY, float nZ,
				 float x1, float y1, float z1,
				 float x2, float y2, float z2,
				 float x3, float y3, float z3
				 ) {
	
	triangles.push_back(STLFace());
	float data[12];
	
	data[0] = nX;
	data[1] = nY;
	data[2] = nZ;
	
	data[3] = x1;
	data[4] = y1;
	data[5] = z1;
	
	data[6] = x2;
	data[7] = y2;
	data[8] = z2;
	
	data[9] = x3;
	data[10] = y3;
	data[11] = z3;
	
	
	triangles[triangles.size()-1].parseFace((char*)data);
	
}

void ofxSTLModel::clear() {
    
    triangles.clear();
    vboMesh.clear();
    
    vertices.clear();
//    normals.clear();
}

