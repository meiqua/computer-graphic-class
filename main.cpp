#ifdef WIN32
	#include <windows.h>
#endif
#include "glut.h"
#include "mesh.h"
#include "matrix.h"
#include "OpenGLProjector.h"
#include <fstream>//20110819
#include "wtypes.h"
#include <map>
// Enumeration
enum EnumDisplayMode { HIDDENLINE, FLATSHADED, SMOOTHSHADED, COLORSMOOTHSHADED,DELETESELECTEDVERTEX,WIREFRAME};

enum Mode 
{ 
	Viewing, 
	Selection,
};
Mode currentMode = Viewing;

// variables
int displayMode = FLATSHADED;	// current display mode
int mainMenu, displayMenu;		// glut menu handlers
int winWidth, winHeight;		// window width and height
double winAspect;				// winWidth / winHeight;
int lastX, lastY;				// last mouse motion position
int currSelectedVertex = -1;         // current selected vertex
bool leftDown,leftUp, rightUp, rightDown, middleDown, middleUp, shiftDown;		// mouse down and shift down flags
double sphi = 90.0, stheta = 45.0, sdepth = 10;	// for simple trackball
double xpan = 0.0, ypan = 0.0, currentR = 10;				// for simple translation
double currentXangle = 0, currentYangle = 0;

double zNear = 1.0, zFar = 100.0;
double g_fov = 45.0;
Vector3d g_center;
double g_sdepth;
Mesh mesh;	// our mesh
int spFaceOnly = 0;

// functions
void SetBoundaryBox(const Vector3d & bmin, const Vector3d & bmax);
void InitGL();
void InitMenu();
void InitGeometry();
void MenuCallback(int value);
void ReshapeFunc(int width, int height);
void DisplayFunc();
void DrawWireframe();
void DrawHiddenLine();
void DrawFlatShaded();
void DrawSmoothShaded();
void DrawColorSmoothShaded();
void DrawSelectedVertices();
//void Partition();

void KeyboardFunc(unsigned char ch, int x, int y);
void MouseFunc(int button, int state, int x, int y);
void MotionFunc(int x, int y);
void SelectVertexByPoint();
void DeleteSelectedVertex(int vertex);

vector<HEdgeList> edgeUnions;
HEdgeList spBoundary;

int numComponents=0;

vector<HEdgeList> uf(HEdgeList & spBoundary) {
	// union find method to find edgeUnions
	vector<HEdgeList> eUnions;
	vector<int> union_find(spBoundary.size());
	for (size_t i = 0; i < union_find.size(); i++) {
		union_find[i] = i;
	}
	for (size_t i = 0; i < spBoundary.size() - 1; i++) {
		int parent = i;
		int child;
		for (size_t i1 = i + 1; i1 < spBoundary.size(); i1++) {
			// if have same vertex, regard as connection
			if (spBoundary[i]->Next()->Start() == spBoundary[i1]->Start() ||
				spBoundary[i]->Next()->Start() == spBoundary[i1]->Next()->Start() ||
				spBoundary[i]->Start() == spBoundary[i1]->Start() ||
				spBoundary[i]->Start() == spBoundary[i1]->Next()->Start()) {
				child = i1;
				vector<int> passBy;
				while (child != union_find[child]) //find child root
				{
					passBy.push_back(child);
					child = union_find[child];
				}
				for (size_t i11 = i + 1; i11 < passBy.size(); i11++) {
					union_find[passBy[i]] = child;
				}
				passBy.clear();
				while (parent != union_find[parent]) //find parent root
				{
					passBy.push_back(parent);
					parent = union_find[parent];
				}
				for (size_t i11 = i + 1; i11 < passBy.size(); i11++) {
					union_find[passBy[i]] = parent;
				}
				union_find[child] = parent;
			}
		}
	}

	vector< vector<int> > unions;
	std::map<int, int> rootToUnions;
	for (size_t i = 0; i < union_find.size(); i++) {
		if (union_find[i] == i) {
			vector<int> aUnion;
			aUnion.push_back(i);
			unions.push_back(aUnion);
			rootToUnions[i] = unions.size() - 1;
		}
	}
	for (size_t i = 0; i < union_find.size(); i++) {
		int parent = union_find[i];
		while (parent != union_find[parent]) //find root parent
		{
			parent = union_find[parent];
		}
		int whichUnion = rootToUnions.find(parent)->second;
		vector<int> aUnion = unions[whichUnion];
		aUnion.push_back(i);
		unions[whichUnion] = aUnion;
	}

	// get edgeUnions
	for (size_t i = 0; i < unions.size(); i++) {
		vector<int> aUnion = unions[i];
		HEdgeList egList;
		for (size_t i1 = 0; i1 < aUnion.size(); i1++) {
			egList.push_back(spBoundary[aUnion[i1]]);
		}
		eUnions.push_back(egList);
	}
	return eUnions;
}





void findSp() {
	FaceList fList = mesh.Faces();
	HEdgeList heList = mesh.Edges();
	HEdgeList bheList = mesh.BoundaryEdges();
	for (size_t i = 0; i<fList.size(); i++) {
		if (fList[i] != NULL && fList[i]->HalfEdge()->LeftFace() != NULL)
		{
			Face *f = fList[i];
			const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
			const Vector3d & pos2 = f->HalfEdge()->End()->Position();
			const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
			Vector3d normal = (pos2 - pos1).Cross(pos3 - pos1);
			normal /= normal.L2Norm();

			f->SetNormal_f(normal);
			if (normal[1] < 0) {
				f->HalfEdge()->SetValid(true);
				f->HalfEdge()->Next()->SetValid(true);
				f->HalfEdge()->Next()->Next()->SetValid(true);
			}
			else {
				f->HalfEdge()->SetValid(false);
				f->HalfEdge()->Next()->SetValid(false);
				f->HalfEdge()->Next()->Next()->SetValid(false);
			}
		}
	}

	for (size_t i = 0; i < heList.size(); i++) {
		HEdge* edge = heList[i];
		if (edge->IsValid()) {
			if ((!edge->Twin()->IsValid())) {
				spBoundary.push_back(edge);
			}
		}
	}
	for (size_t i = 0; i < bheList.size(); i++) {
		HEdge* edge = bheList[i];
		if (edge->IsValid()) {
			spBoundary.push_back(edge);
		}
	}
//	edgeUnions = uf(spBoundary);
}


void SetBoundaryBox(const Vector3d & bmin, const Vector3d & bmax) {
	double PI = 3.14159265358979323846;
	double radius = bmax.Distance(bmin);
	g_center = 0.5 * (bmin+bmax);
	zNear    = 0.2 * radius / sin(0.5 * g_fov * PI / 180.0);
	zFar     = zNear + 2.0 * radius;
	g_sdepth = zNear + radius;
	zNear *= 0.1;
	zFar *= 10;
	sdepth = g_sdepth;
	currentR = sdepth;
}

// init openGL environment
void InitGL() {
	GLfloat light0Position[] = { 0, 1, 0, 1.0 }; 

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutCreateWindow("Mesh Viewer");
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glPolygonOffset(1.0, 1.0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glLightfv (GL_LIGHT0, GL_POSITION, light0Position);
	glEnable(GL_LIGHT0);

	glutReshapeFunc(ReshapeFunc);
	glutDisplayFunc(DisplayFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
}

// init right-click menu
void InitMenu() {
	displayMenu = glutCreateMenu(MenuCallback);
	glutAddMenuEntry("Wireframe", WIREFRAME);
	glutAddMenuEntry("Hidden Line", HIDDENLINE);
	glutAddMenuEntry("Flat Shaded", FLATSHADED);
	glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);
	glutAddMenuEntry("Color Smooth Shaded", COLORSMOOTHSHADED);
    glutAddMenuEntry("Delete Selected Vertex", DELETESELECTEDVERTEX);
	mainMenu = glutCreateMenu(MenuCallback);
	glutAddSubMenu("Display", displayMenu);
	glutAddMenuEntry("Exit", 99);
	glutAttachMenu(GLUT_RIGHT_BUTTON);//glutAttachMenu(GLUT_MIDDLE_BUTTON);
}

// init geometry (if no input argument is provided)
void InitGeometry() {
	const int VSIZE = 4;
	const int HESIZE = 12;
	const int FSIZE = 4;
	int i;
	Vertex *v[VSIZE];
	HEdge *he[HESIZE];
	Face *f[FSIZE];
	
	for (i=0; i<VSIZE; i++) {
		v[i] = new Vertex();
		mesh.vList.push_back(v[i]);
	}
	v[0]->SetPosition(Vector3d(0.0, 0.0, 0.0));
	v[1]->SetPosition(Vector3d(10.0, 0.0, 0.0));
	v[2]->SetPosition(Vector3d(0.0, 10.0, 0.0));
	v[3]->SetPosition(Vector3d(0.0, 0.0, 10.0));

	v[0]->SetNormal(Vector3d(-0.577, -0.577, -0.577));
	v[1]->SetNormal(Vector3d(0.0, -0.7, -0.7));
	v[2]->SetNormal(Vector3d(-0.7, 0.0, -0.7));
	v[3]->SetNormal(Vector3d(-0.7, -0.7, 0.0));

	for (i=0; i<FSIZE; i++) {
		f[i] = new Face();
		mesh.fList.push_back(f[i]);
	}

	for (i=0; i<HESIZE; i++) {
		he[i] = new HEdge();
		mesh.heList.push_back(he[i]);
	}
	for (i=0; i<FSIZE; i++) {
		int base = i*3;
		SetPrevNext(he[base], he[base+1]);
		SetPrevNext(he[base+1], he[base+2]);
		SetPrevNext(he[base+2], he[base]);
		SetFace(f[i], he[base]);
	}
	SetTwin(he[0], he[4]);
	SetTwin(he[1], he[7]);
	SetTwin(he[2], he[10]);
	SetTwin(he[3], he[8]);
	SetTwin(he[5], he[9]);
	SetTwin(he[6], he[11]);
	he[0]->SetStart(v[1]); he[1]->SetStart(v[2]); he[2]->SetStart(v[3]);
	he[3]->SetStart(v[0]); he[4]->SetStart(v[2]); he[5]->SetStart(v[1]);
	he[6]->SetStart(v[0]); he[7]->SetStart(v[3]); he[8]->SetStart(v[2]);
	he[9]->SetStart(v[0]); he[10]->SetStart(v[1]); he[11]->SetStart(v[3]);
	v[0]->SetHalfEdge(he[3]);
	v[1]->SetHalfEdge(he[0]);
	v[2]->SetHalfEdge(he[1]);
	v[3]->SetHalfEdge(he[2]);
}

// GLUT menu callback function
void MenuCallback(int value) {
	switch (value) {
	case 99: exit(0); break;
	default: 
		displayMode = value;
		glutPostRedisplay();
		break;
	}
}

// GLUT reshape callback function
void ReshapeFunc(int width, int height) {
	winWidth = width;
	winHeight = height;
	winAspect = (double)width/(double)height;
	glViewport(0, 0, width, height);
	glutPostRedisplay();
}

// GLUT display callback function
void DisplayFunc() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(g_fov, winAspect, zNear, zFar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glRotatef(currentYangle, 1.0, 0.0, 0.0);
	glRotatef(-currentXangle, 0, 1.0, 0.0);
	glTranslatef(0, 0, -currentR/cos(currentYangle/180*3.14)/ cos(currentXangle/180*3.14));

	glRotatef(-stheta, 1.0, 0.0, 0.0);
	glRotatef(sphi, 0.0, 1.0, 0.0);
	glTranslatef(-g_center[0], -g_center[1], -g_center[2]);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	switch (displayMode) {
	case WIREFRAME: DrawWireframe(); break;
	case HIDDENLINE: DrawHiddenLine(); break;
	case FLATSHADED: DrawFlatShaded(); break;
	case SMOOTHSHADED: DrawSmoothShaded(); break;
	case COLORSMOOTHSHADED: DrawColorSmoothShaded(); break;
	case DELETESELECTEDVERTEX: DeleteSelectedVertex(currSelectedVertex); break;
	}

	glutSwapBuffers();
}

// Wireframe render function
void DrawWireframe() {
	HEdgeList heList = mesh.Edges();
	HEdgeList bheList = mesh.BoundaryEdges();
	glColor3f(0.3, 0.3, 1.0);
	glBegin(GL_LINES);
	size_t i;
	for (i=0; i<heList.size(); i++) {
		glVertex3dv(heList[i]->Start()->Position().ToArray());
		glVertex3dv(heList[i]->End()->Position().ToArray());
	}
    
	glColor3f(1, 0, 0);
	for (i=0; i<bheList.size(); i++) {
		glVertex3dv(bheList[i]->Start()->Position().ToArray());
		glVertex3dv(bheList[i]->End()->Position().ToArray());
	}
	
	glEnd();
}

// Hidden Line render function
void DrawHiddenLine() {
	FaceList fList = mesh.Faces();
	glShadeModel(GL_FLAT); 
	glEnable(GL_POLYGON_OFFSET_FILL);
	glColor3f(0, 0, 0);
	glBegin(GL_TRIANGLES);
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
		const Vector3d & pos2 = f->HalfEdge()->End()->Position();
		const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
		glVertex3dv(pos1.ToArray());
		glVertex3dv(pos2.ToArray());
		glVertex3dv(pos3.ToArray());
	}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);

	DrawWireframe();
}

// Flat Shaded render function
void DrawFlatShaded() {
	FaceList fList = mesh.Faces();
	glShadeModel(GL_FLAT); 
	glEnable(GL_LIGHTING);
	glColor3f(0.4f, 0.4f, 1.0f);
	glBegin(GL_TRIANGLES);
	for (size_t i=0; i<fList.size(); i++) {
		if(fList[i]!=NULL && fList[i]->HalfEdge()->LeftFace()!=NULL )
		{
		Face *f = fList[i];
		const Vector3d & pos1 = f->HalfEdge()->Start()->Position();
		const Vector3d & pos2 = f->HalfEdge()->End()->Position();
		const Vector3d & pos3 = f->HalfEdge()->Next()->End()->Position();
		Vector3d normal = (pos2-pos1).Cross(pos3-pos1);
		normal /= normal.L2Norm();
        
		f->SetNormal_f(normal);//1007
        
		if (spFaceOnly%3==0) {
			glNormal3dv(normal.ToArray());
			glVertex3dv(pos1.ToArray());
			glVertex3dv(pos2.ToArray());
			glVertex3dv(pos3.ToArray());
		}
		else if (spFaceOnly % 3 == 1) {
			if (normal[1] < 0) {
				glNormal3dv(normal.ToArray());
				glVertex3dv(pos1.ToArray());
				glVertex3dv(pos2.ToArray());
				glVertex3dv(pos3.ToArray());
			}
		}
		}
	}
	glEnd();
	glDisable(GL_LIGHTING);

	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	size_t i;
	for (i = 0; i<spBoundary.size(); i++) {
		glVertex3dv(spBoundary[i]->Start()->Position().ToArray());
		glVertex3dv(spBoundary[i]->End()->Position().ToArray());
	}
	glEnd();
}

// Smooth Shaded render function
void DrawSmoothShaded() { 
	FaceList fList = mesh.Faces();
	glShadeModel(GL_SMOOTH); 
	glEnable(GL_LIGHTING);
	glColor3f(0.4f, 0.4f, 1.0f);
	glBegin(GL_TRIANGLES) ;
	for (size_t i=0; i<fList.size(); i++) {
		Face *f = fList[i];
		Vertex * v1 = f->HalfEdge()->Start();
		Vertex * v2 = f->HalfEdge()->End();
		Vertex * v3 = f->HalfEdge()->Next()->End();
		glNormal3dv(v1->Normal().ToArray());
		glVertex3dv(v1->Position().ToArray());
		glNormal3dv(v2->Normal().ToArray());
		glVertex3dv(v2->Position().ToArray());
		glNormal3dv(v3->Normal().ToArray());
		glVertex3dv(v3->Position().ToArray());
	}
	glEnd();
	glDisable(GL_LIGHTING);
}

void DrawColorSmoothShaded() {
	cout<< "the colored smooth model is drawn"<<endl;
}

 

// draw the selected ROI vertices on the mesh
void DrawSelectedVertices()
{
cout<< "the selected vertex is drawn "<<endl;
}

//delete selected vertex and its incident faces and half-edges

void DeleteSelectedVertex(int vertex)
{
   
	cout<< "the selected vertex is gone, its index is "<<endl;
} 


// GLUT keyboard callback function
void KeyboardFunc(unsigned char ch, int x, int y) { 
	switch (ch) { 
	case '3':
	
		cout<<"Translation mode"<<endl;

		break; 
	case '4':
			spFaceOnly++;
		cout<<"0004"<<endl;;

		break;
	case '5':
        cout<<"0005"<<endl;
        break;
	case '6':
        cout<<"0006"<<endl;
        break;

	case '1':	// key '1'
		currentMode = Viewing;
		cout << "Viewing mode" << endl;
		break;
	case '2':	// key '2'
		currentMode = Selection;
		cout << "Selection mode" << endl;
		break;

	case '9': 
		cout<<"0009"<<endl;
		break;
		//cout<<"test000 "<<currSelectedVertex<<endl;
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

// GLUT mouse callback function
void MouseFunc(int button, int state, int x, int y) {
	
	lastX = x;
	lastY = y;
	leftDown = (button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN);
	leftUp = (button == GLUT_LEFT_BUTTON) && (state == GLUT_UP);
    rightDown = (button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN);
	rightUp = (button == GLUT_RIGHT_BUTTON) && (state == GLUT_UP);
	middleDown = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN);
	middleUp = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_UP);
	shiftDown = (glutGetModifiers() & GLUT_ACTIVE_SHIFT);

	if (currentMode == Selection && state == GLUT_UP)
	{
	    if (middleUp)
		{
		    if (currSelectedVertex != -1)
		    {
		        mesh.Vertices()[currSelectedVertex]->SetFlag(0);
		        currSelectedVertex = -1;
		    }
		}
		else SelectVertexByPoint();
		
		//if (leftUp)
	    //{
		  //  SelectVertexByPoint(); 
			
		//}
        //if(middleUp)
		//{DeleteSelectedVertex(currSelectedVertex);}
        
		lastX = lastY = 0;
		glutPostRedisplay();
	}

}

// GLUT mouse motion callback function
void MotionFunc(int x, int y) {
	if (leftDown)
		if(!shiftDown) { // rotate
			sphi += (double)(x - lastX) / 4.0;
			stheta += (double)(lastY - y) / 4.0;
		} else { // translate
//********* please input your tranlation code here *************//
			currentXangle += double(x - lastX) / zNear / 10000 /3.14*180;
			currentYangle += double(lastY - y) / zNear / 10000 /3.14*180;
//**************************************************************//
		}
	// scale
	if (middleDown) currentR += (double)(lastY - y) / 10.0;

	lastX = x;
	lastY = y;
	glutPostRedisplay();
}


// select a mesh point
void SelectVertexByPoint()
{
	// get the selection position
	int x = lastX, y = winHeight - lastY;
	Vector3d u(x,y,0);
	
	OpenGLProjector projector; 

	
}



// main function
void main(int argc, char **argv) {
	glutInit(&argc, argv);
	InitGL();
	InitMenu();
	if (argc>=2) mesh.LoadObjFile(argv[1]);
	else InitGeometry();
	SetBoundaryBox(mesh.MinCoord(), mesh.MaxCoord());
	/************************************************************************/
	/* activate the following code if you finish the corresponding functions*/
 	mesh.DisplayMeshInfo();
	/************************************************************************/
	findSp();
	glutMainLoop();
}


