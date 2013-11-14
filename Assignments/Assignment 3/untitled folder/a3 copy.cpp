#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <map>

#define LEFT_EYE 0
#define RIGHT_EYE 1
#define NOSE 2
#define LEFT_MOUTH 3
#define RIGHT_MOUTH 4
#define CHIN 5

#define PART_COUNT 6
#define STEP_SIZE 7

using namespace std;

// Definition for a 2D pixel point.
// Add additional fields if necessary.
// 
class Coordinate
{  
public:
  Coordinate() :row(0), col(0) {};
  double row, col;
};

// Ground truth nformation for training or testing image
class GTImage
{
public:
  string image_filename;
  Coordinate coords[6];    
};

// Code for computing unary potential functions
//  Feel free to modify this if you want, but 
//  it should work as-is.
#include <UnaryPotential.h>


// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlay rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}





// Load ground truth coordinates and image names from file
//
vector<GTImage> load_groundtruth(const string &gt_filename)
{
  // read in training images, one at a time, and put the ground truth coordinates
  //  into a vector.
  ifstream ifs(gt_filename.c_str());
  
  vector<GTImage> training_images;
  while(ifs.good())
    {
      string image_number;
      ifs >> image_number;
      if(!ifs.good())	break;
      
      GTImage gt;
      gt.image_filename = string("imgs/image_") + image_number + ".png";
      for(int ii=0; ii<6; ii++)
	ifs >> gt.coords[ii].col >> gt.coords[ii].row;

      training_images.push_back(gt);
    }

  return training_images;
}



// Draw a given part configuration on an image, using simple
//  colored rectangles.
// Color code is: 
// left eye: bright green
// right eye: dark green
// nose: bright red
// left mouth: bright blue
// right mouth: dark blue
// chin: dark red
void draw_configuration(const char *output_fname, const char *sample_filename, const vector<Coordinate> &configuration)
{
  int g_code[] = {255, 128,   0,   0,   0, 0};
  int r_code[] = {0, 0,     255,   0,   0, 128};
  int b_code[] = {0, 0,       0, 255, 128, 0};

  SDoublePlane img = SImageIO::read_png_file(sample_filename);
  SDoublePlane r = img, g = img, b = img;
  for(int ii=0; ii<configuration.size(); ii++)
    {
      overlay_rectangle(g, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, g_code[ii], 3);
      overlay_rectangle(r, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, r_code[ii], 3);
      overlay_rectangle(b, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, b_code[ii], 3);
    }

  SImageIO::write_png_file(output_fname, r, g, b);
}

/*
The following function calculates the mean and variance of (Li-Lj) for all images in the training data set, where Li,Lj are the parts in the image
*/
void learning(vector<GTImage> training_data, Coordinate *means[6], Coordinate *vars[6])
{
  
    //intialize all values to zero
    for (int k=0; k<PART_COUNT; k++)
    {
        for (int j=k+1; j<PART_COUNT; j++) //means for combinations of parts, eg Mu: righteye,lefteye would be means[0][1]
        {                                // Mu used in message right eye to left will contain negated version of mu in [0][1]               
          means[k][j].row=0;                  //and be store in [1][0]
          means[k][j].col=0;
          vars[k][j].row=0;
          vars[k][j].col=0;            //variances set up same way as means
        }
    }
    //sum up differences of combinations of parts to calculate means for them
    for (int i=0; i<training_data.size(); i++)
    {
      for (int k=0; k<PART_COUNT; k++)
      {
        for (int j=k+1; j<PART_COUNT; j++)
        {
          means[k][j].row += (training_data[i].coords[k].row - training_data[i].coords[j].row) ;
          means[k][j].col += (training_data[i].coords[k].col - training_data[i].coords[j].col) ;       
        }
      }
    }
    //caclulate mean of differences by dividing sums by total number of elements
    for (int k=0; k<PART_COUNT; k++)
    {
    	for (int j=k+1; j<PART_COUNT; j++)
    	{
    		means[k][j].row = means[k][j].row/training_data.size();
    		means[k][j].col = means[k][j].col/training_data.size();
    	}
    }
    
    //calculations of sum of (observations - corresponding mean)sqaured to calculate variances 
    for (int i=0; i<training_data.size(); i++)
    {
    	for (int k=0; k<PART_COUNT; k++)
    	{
    		for(int j=k+1; j<PART_COUNT; j++)
    		{

          vars[k][j].row +=  pow( ( (training_data[i].coords[k].row - training_data[i].coords[j].row) - means[k][j].row) , 2);
          vars[k][j].col +=  pow( ( (training_data[i].coords[k].col- training_data[i].coords[j].col) - means[k][j].col) , 2);
    		}
    	}
    }
  
  //divide sums to get variances 
  for (int k=0; k<PART_COUNT; k++)
  {
  	for (int j=k+1; j<PART_COUNT; j++)
  	{
  		vars[k][j].row = vars[k][j].row/training_data.size();
  		vars[k][j].col = vars[k][j].col/training_data.size();
      //cout << "Means of " << k << "  " << j << " " << means[k][j].row << " " <<  means[k][j].col << endl;
      //cout << "Vars of " << k << "  " << j << " " << vars[k][j].row << " " <<  vars[k][j].col << endl;
  	}
  }
  
}


//its called Phi_calculator but its actually calculating deltas (messages) since we are adding the unary of the child in question
SDoublePlane Phi_calculator(vector <vector<Coordinate> > &coor_map, Coordinate mean, Coordinate var, SDoublePlane my_potentials) //returns an array of coordinates, which is the best possible assignments for two parts i and j
{
  int rows = my_potentials.rows();
  int cols = my_potentials.cols();
  cout << "eac here " << endl;
  //cout << "rows " << rows <<  " cols " << cols << endl;
  SDoublePlane phi = SDoublePlane(rows,cols);  //initialize potential table to really low value
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      phi[i][j] = -1e100;
    }

  }
  for(int ix = 0; ix<rows; ix+=STEP_SIZE)    //HEART OF THE COMPUTATION DONE BY THIS PROGRAM
  {
    for(int iy = 0; iy<cols; iy+=STEP_SIZE) //outer two loops iterate all possible coordinates for part i
    {
      double max_phi = -1e100;
      Coordinate max_coord;
      max_coord.row = 0;
      max_coord.col = 0;
      for (int jx = 0; jx<rows; jx+=STEP_SIZE)      //inner two loops iterate through all possible coordinates of part j, for each of part i's value
      {
        for (int jy = 0; jy<cols; jy+=STEP_SIZE)
        {
          double t1 = pow(ix-jx-mean.row,2)/(2*var.row);
          double t2 = pow(iy-jy-mean.col,2)/(2*var.col);    //two parts of the gaussian function from assignment document
          double temp = -t1-t2 + my_potentials[jx][jy];     //add unarypotential of part j. This turns our calculation into a delta function. 
          if(temp>max_phi)                
          {
            max_phi = temp;
            max_coord.row = jx;       //maximize and store for each possible assignment of parent (which is iterating in outer loop)
            max_coord.col = jy;       //a coordinate value of part j which is most likely assignment for that value of part i
                                      //as well as a delta value 
          }
        }
      }
      phi[ix][iy] = max_phi;  //the two tables phi and coor_map are essentially mappings of potential values to coordinates
      coor_map[ix][iy].row = max_coord.row;  //where in the phi DoublePlane index [4][15] we store a potential value that is the highest potential value                                        //for part i's being at 4,15 out of part j's all possible values
      coor_map[ix][iy].col = max_coord.col;  //similar coor_map index 4,15 we store for the part i's being at 4,15, what is the part j's most likely coordinate
                                             //assignment
    }

  }

  return phi; //this double plane storing potentials is returned. 
             //Although we want to return two tables, the other one (coor_map) is being passed by reference
            //and its values are simply updated.

} 





//method for performing MAP inference on model B, takes potentials
GTImage MAP_Inference_B (vector<SDoublePlane> my_potentials, string filename, Coordinate *means[6], Coordinate *vars[6]) 
{
//takes all training data loaded by load_ground_truth and means and vars meade by learning(). Also takes mypotentials which are already ready in main from where this method will be called.

  int rows = my_potentials[0].rows();//number of rows in the image
  int cols = my_potentials[0].cols();	//number of cols in the image
  
  //********************************************GRAPH B *********************************************************
  vector < vector <Coordinate> > coor_map_le_re(rows,vector<Coordinate>(cols)); // declare these vectors of vectors of coordinates
  vector < vector <Coordinate> > coor_map_le_n(rows,vector<Coordinate>(cols));  // to store most likely coordinate of child in the index
  vector < vector <Coordinate> > coor_map_le_rm(rows,vector<Coordinate>(cols)); // which denotes coordinate of parent
  vector < vector <Coordinate> > coor_map_le_lm(rows,vector<Coordinate>(cols)); // which in all 5 cases for this model is left eye
  vector < vector <Coordinate> > coor_map_le_c(rows,vector<Coordinate>(cols));
  for (int ii = 0;ii<rows;ii++)
  {
    for(int jj = 0;jj<cols;jj++)     //intialize these vectors to zero
    {
      coor_map_le_re[ii][jj].row =0;     
      coor_map_le_re[ii][jj].col = 0;

      coor_map_le_n[ii][jj].row =0;
      coor_map_le_n[ii][jj].col = 0;

      coor_map_le_rm[ii][jj].row =0;
      coor_map_le_rm[ii][jj].col = 0;

      coor_map_le_lm[ii][jj].row =0;
      coor_map_le_lm[ii][jj].col = 0;

      coor_map_le_c[ii][jj].row =0;
      coor_map_le_c[ii][jj].col = 0;
    }
  }
  
  //call Phi_calculator for each clique which returns doubleplanes, and updates 2d vectors declared above
  SDoublePlane lefteye_righteye = Phi_calculator(coor_map_le_re, means[LEFT_EYE][RIGHT_EYE], vars[LEFT_EYE][RIGHT_EYE], my_potentials[RIGHT_EYE]);
  SDoublePlane lefteye_nose = Phi_calculator(coor_map_le_n, means[LEFT_EYE][NOSE], vars[LEFT_EYE][NOSE], my_potentials[NOSE]);
  SDoublePlane lefteye_leftmouth = Phi_calculator(coor_map_le_lm, means[LEFT_EYE][LEFT_MOUTH], vars[LEFT_EYE][LEFT_MOUTH], my_potentials[LEFT_MOUTH]);
  SDoublePlane lefteye_rightmouth = Phi_calculator(coor_map_le_rm, means[LEFT_EYE][RIGHT_MOUTH], vars[LEFT_EYE][RIGHT_MOUTH], my_potentials[RIGHT_MOUTH]);
  SDoublePlane lefteye_chin = Phi_calculator(coor_map_le_c, means[LEFT_EYE][CHIN], vars[LEFT_EYE][CHIN], my_potentials[CHIN]);
 
  GTImage answer; //declare instance of GT image which will be return value
  answer.image_filename = filename; 
  vector<Coordinate> assignments_B(6);  //declare vector in case we want to call draw_configuration
  double temp = 0.0;
  double max = -1e100;
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++) //simply iterate through row and columsn, adding up unary of root(left eye) + deltas of its children and maximize
    {
      temp = lefteye_righteye[i][j] + lefteye_nose[i][j] + lefteye_leftmouth[i][j] + lefteye_rightmouth[i][j] + lefteye_chin[i][j] + my_potentials[0][i][j];
      if (temp > max)
      {
        max = temp;
        assignments_B[LEFT_EYE].row = i;
        assignments_B[LEFT_EYE].col = j;
        answer.coords[LEFT_EYE].row = i;
        answer.coords[LEFT_EYE].col = j;  //index at which we maximize root is correct assignment for root
        //assignments_D and answer.coords are both same thing, one is used for draw config
      }
    }
  } 


  assignments_B[RIGHT_EYE].row = coor_map_le_re[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].row;
  assignments_B[RIGHT_EYE].col = coor_map_le_re[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].col;
 
  assignments_B[NOSE].row = coor_map_le_n[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].row;
  assignments_B[NOSE].col = coor_map_le_n[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].col;

  assignments_B[LEFT_MOUTH].row = coor_map_le_lm[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].row;
  assignments_B[LEFT_MOUTH].col = coor_map_le_lm[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].col;

  assignments_B[RIGHT_MOUTH].row = coor_map_le_rm[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].row;
  assignments_B[RIGHT_MOUTH].col = coor_map_le_rm[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].col;

  assignments_B[CHIN].row = coor_map_le_c[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].row;
  assignments_B[CHIN].col = coor_map_le_c[assignments_B[LEFT_EYE].row][assignments_B[LEFT_EYE].col].col;
  
  //assignments to coords of GTImage that will be returned. As you can see they are acquired by using assignment of root (left eye) as index and
  //using it on appropriate coor_map
  answer.coords[RIGHT_EYE].row = coor_map_le_re[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].row;
  answer.coords[RIGHT_EYE].col = coor_map_le_re[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].col;
 
  answer.coords[NOSE].row = coor_map_le_n[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].row;
  answer.coords[NOSE].col = coor_map_le_n[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].col;

  answer.coords[LEFT_MOUTH].row = coor_map_le_lm[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].row;
  answer.coords[LEFT_MOUTH].col = coor_map_le_lm[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].col;

  answer.coords[RIGHT_MOUTH].row = coor_map_le_rm[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].row;
  answer.coords[RIGHT_MOUTH].col = coor_map_le_rm[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].col;

  answer.coords[CHIN].row = coor_map_le_c[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].row;
  answer.coords[CHIN].col = coor_map_le_c[answer.coords[LEFT_EYE].row][answer.coords[LEFT_EYE].col].col;

  int strtIndx=filename.find("_"), endIndx=filename.find(".");
  string imageNumber=filename.substr(strtIndx+1,endIndx-strtIndx-1);
  string imgOutput="B_Output_" + imageNumber + ".png";
  draw_configuration(imgOutput.c_str(), filename.c_str(), assignments_B);
  return answer;
}


//method for performing MAP inference on model C, takes potentials
GTImage MAP_Inference_C (vector<SDoublePlane> my_potentials, string filename, Coordinate *means[6], Coordinate *vars[6]) 
{

  
  int rows = my_potentials[0].rows();//number of rows in the image
  int cols = my_potentials[0].cols(); //number of cols in the image
  //************************************* GRAPH C ***************************************************
  vector < vector <Coordinate> > coor_map_n_c(rows,vector<Coordinate>(cols));  // declare these vectors of vectors of coordinates
  vector < vector <Coordinate> > coor_map_n_rm(rows,vector<Coordinate>(cols)); // to store most likely coordinate of child in the index    
  vector < vector <Coordinate> > coor_map_n_re(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_rm_lm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_re_le(rows,vector<Coordinate>(cols));
  /*   
          nose
        /  |   \
    reye chin rmouth
      |        |
    leye     lmouth     */
      
  for (int ii = 0;ii<rows;ii++)
  {
    for(int jj = 0;jj<cols;jj++)
    {
      coor_map_re_le[ii][jj].row =0;
      coor_map_re_le[ii][jj].col = 0;

      coor_map_rm_lm[ii][jj].row =0;
      coor_map_rm_lm[ii][jj].col = 0;

      coor_map_n_re[ii][jj].row =0;
      coor_map_n_re[ii][jj].col = 0;

      coor_map_n_rm[ii][jj].row =0;
      coor_map_n_rm[ii][jj].col = 0;

      coor_map_n_c[ii][jj].row =0;
      coor_map_n_c[ii][jj].col = 0;
    }
  }

  means[NOSE][RIGHT_EYE].row=-means[RIGHT_EYE][NOSE].row;
  means[NOSE][RIGHT_EYE].col=-means[RIGHT_EYE][NOSE].col;
  
  vars[NOSE][RIGHT_EYE].row=vars[RIGHT_EYE][NOSE].row;
  vars[NOSE][RIGHT_EYE].col=vars[RIGHT_EYE][NOSE].col;
  //These double planes can be thought of as messages
  //after calls in two lines below, these two double planes will contain delta of pairwise potentials of rm_lm + unary of lm 

  SDoublePlane rightmouth_leftmouth = Phi_calculator(coor_map_rm_lm, means[LEFT_MOUTH][RIGHT_MOUTH], vars[LEFT_MOUTH][RIGHT_MOUTH], my_potentials[LEFT_MOUTH]); //pairwise of lm,rm + unary lm
  SDoublePlane righteye_lefteye = Phi_calculator(coor_map_re_le, means[LEFT_EYE][RIGHT_EYE], vars[LEFT_EYE][RIGHT_EYE], my_potentials[LEFT_EYE]);  //pairwise of re,le + unary le

  //now we add unary's of rmouth and reye to the two SDoublePLanes above 
  for (int i=0; i< rows; i++)
  {
    for (int j=0; j< cols; j++)
    {
      rightmouth_leftmouth[i][j] += my_potentials[RIGHT_MOUTH][i][j]; //add unary of rm
      righteye_lefteye[i][j] += my_potentials[RIGHT_EYE][i][j]; //add unary re
    }
  } 
  //now right eye and right mouth want to send a message to nose, but instead of unary potentials as parameters they will send two doubleplanes above
  //as paraments since right eye and right mouth themselves have children
  SDoublePlane nose_righteye = Phi_calculator(coor_map_n_re, means[NOSE][RIGHT_EYE], vars[NOSE][RIGHT_EYE], righteye_lefteye);
  SDoublePlane nose_rightmouth = Phi_calculator(coor_map_n_rm, means[NOSE][RIGHT_MOUTH], vars[NOSE][RIGHT_MOUTH], rightmouth_leftmouth);
  //chin has no children so it will only send its unary potentials 
  SDoublePlane nose_chin = Phi_calculator(coor_map_n_c, means[NOSE][CHIN], vars[NOSE][CHIN], my_potentials[CHIN]);


  GTImage answer; //declare instance of GT image 
  answer.image_filename = filename; 
  vector<Coordinate> assignments_C(6);

  double temp = 0.0;
  double max = -1e100;

  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      temp = nose_rightmouth[i][j] + nose_righteye[i][j] + nose_chin[i][j] + my_potentials[2][i][j];
      if (temp > max)
      {
        max = temp;
        assignments_C[NOSE].row = i;
        assignments_C[NOSE].col = j;
        answer.coords[NOSE].row = i;
        answer.coords[NOSE].col = j;      //here we assign values to root
        //assignments_D and answer.coords are both same thing, one is used for draw config
      }
    }
  } 

  //assignments to coords of GTImage that will be returned. As you can see they are acquired by using assignment of root (left eye) as index and
  //using it on appropriate coor_map
  answer.coords[RIGHT_EYE].row = coor_map_n_re[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[RIGHT_EYE].col = coor_map_n_re[answer.coords[NOSE].row][answer.coords[NOSE].col].col;
 
  answer.coords[RIGHT_MOUTH].row = coor_map_n_rm[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[RIGHT_MOUTH].col = coor_map_n_rm[answer.coords[NOSE].row][answer.coords[NOSE].col].col;

  answer.coords[CHIN].row = coor_map_n_c[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[CHIN].col = coor_map_n_c[answer.coords[NOSE].row][answer.coords[NOSE].col].col;

  answer.coords[LEFT_EYE].row = coor_map_re_le[answer.coords[RIGHT_EYE].row][answer.coords[RIGHT_EYE].col].row;
  answer.coords[LEFT_EYE].col = coor_map_re_le[answer.coords[RIGHT_EYE].row][answer.coords[RIGHT_EYE].col].col;

  answer.coords[LEFT_MOUTH].row = coor_map_rm_lm[answer.coords[RIGHT_MOUTH].row][answer.coords[RIGHT_MOUTH].col].row;
  answer.coords[LEFT_MOUTH].col = coor_map_rm_lm[answer.coords[RIGHT_MOUTH].row][answer.coords[RIGHT_MOUTH].col].col;

  assignments_C[RIGHT_EYE].row = coor_map_n_re[assignments_C[NOSE].row][assignments_C[NOSE].col].row;
  assignments_C[RIGHT_EYE].col = coor_map_n_re[assignments_C[NOSE].row][assignments_C[NOSE].col].col;
 
  assignments_C[RIGHT_MOUTH].row = coor_map_n_rm[assignments_C[NOSE].row][assignments_C[NOSE].col].row;
  assignments_C[RIGHT_MOUTH].col = coor_map_n_rm[assignments_C[NOSE].row][assignments_C[NOSE].col].col;

  assignments_C[CHIN].row = coor_map_n_c[assignments_C[NOSE].row][assignments_C[NOSE].col].row;
  assignments_C[CHIN].col = coor_map_n_c[assignments_C[NOSE].row][assignments_C[NOSE].col].col;

  assignments_C[LEFT_EYE].row = coor_map_re_le[assignments_C[RIGHT_EYE].row][assignments_C[RIGHT_EYE].col].row;
  assignments_C[LEFT_EYE].col = coor_map_re_le[assignments_C[RIGHT_EYE].row][assignments_C[RIGHT_EYE].col].col;

  assignments_C[LEFT_MOUTH].row = coor_map_rm_lm[assignments_C[RIGHT_MOUTH].row][assignments_C[RIGHT_MOUTH].col].row;
  assignments_C[LEFT_MOUTH].col = coor_map_rm_lm[assignments_C[RIGHT_MOUTH].row][assignments_C[RIGHT_MOUTH].col].col;

  
  int strtIndx=filename.find("_"), endIndx=filename.find(".");
  string imageNumber=filename.substr(strtIndx+1,endIndx-strtIndx-1);
  string imgOutput="C_Output_" + imageNumber + ".png";
  draw_configuration(imgOutput.c_str(), filename.c_str(), assignments_C);
  
  return answer; //returns answers in a nice GTImage format to make it easy to compare in evaluation

} 
GTImage MAP_Inference_D (vector<SDoublePlane> my_potentials, string filename, Coordinate *means[6], Coordinate *vars[6]) 
{
  
  
  int rows = my_potentials[0].rows();//number of rows in the image
  int cols = my_potentials[0].cols(); //number of cols in the image
  //***********************************************GRAPH D****************************************************
  vector < vector <Coordinate> > coor_map_n_le(rows,vector<Coordinate>(cols)); // declare these vectors of vectors of coordinates
  vector < vector <Coordinate> > coor_map_n_re(rows,vector<Coordinate>(cols)); // to store most likely coordinate of child in the index  
  vector < vector <Coordinate> > coor_map_n_lm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_rm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_c(rows,vector<Coordinate>(cols));
  //initialize 
  for (int ii = 0;ii<rows;ii++) //initialize coordinate maps to zero
  {
    for(int jj = 0;jj<cols;jj++)
    {
      coor_map_n_le[ii][jj].row =0;
      coor_map_n_le[ii][jj].col = 0;

      coor_map_n_re[ii][jj].row =0;
      coor_map_n_re[ii][jj].col = 0;

      coor_map_n_lm[ii][jj].row =0;
      coor_map_n_lm[ii][jj].col = 0;

      coor_map_n_rm[ii][jj].row =0;
      coor_map_n_rm[ii][jj].col = 0;

      coor_map_n_c[ii][jj].row =0;
      coor_map_n_c[ii][jj].col = 0;
    }
  }
  //These double planes can be thought of as messages. We negate them because message is going from nose to each of the eyes
  means[NOSE][LEFT_EYE].row = -means[LEFT_EYE][NOSE].row;
  means[NOSE][LEFT_EYE].col = -means[LEFT_EYE][NOSE].col;
  means[NOSE][RIGHT_EYE].row = -means[RIGHT_EYE][NOSE].row;
  means[NOSE][RIGHT_EYE].col = -means[RIGHT_EYE][NOSE].col;

  vars[NOSE][LEFT_EYE].row = vars[LEFT_EYE][NOSE].row;
  vars[NOSE][LEFT_EYE].col = vars[LEFT_EYE][NOSE].col;
  vars[NOSE][RIGHT_EYE].row = vars[RIGHT_EYE][NOSE].row;
  vars[NOSE][RIGHT_EYE].col = vars[RIGHT_EYE][NOSE].col;
  
  SDoublePlane nose_lefteye = Phi_calculator(coor_map_n_le, means[NOSE][LEFT_EYE], vars[NOSE][LEFT_EYE], my_potentials[LEFT_EYE]);
  SDoublePlane nose_righteye = Phi_calculator(coor_map_n_re, means[NOSE][RIGHT_EYE], vars[NOSE][RIGHT_EYE], my_potentials[RIGHT_EYE]);
  SDoublePlane nose_leftmouth = Phi_calculator(coor_map_n_lm, means[NOSE][LEFT_MOUTH], vars[NOSE][LEFT_MOUTH], my_potentials[LEFT_MOUTH]);
  SDoublePlane nose_rightmouth = Phi_calculator(coor_map_n_rm, means[NOSE][RIGHT_MOUTH], vars[NOSE][RIGHT_MOUTH], my_potentials[RIGHT_MOUTH]);
  SDoublePlane nose_chin = Phi_calculator(coor_map_n_c, means[NOSE][CHIN], vars[NOSE][CHIN], my_potentials[CHIN]);

  GTImage answer; //declare instance of GT image which will be return value
  answer.image_filename = filename; 
  SDoublePlane nose = SDoublePlane(rows, cols);
  vector<Coordinate> assignments_D(6); 
  double temp = 0.0;
  double max = -1e100;
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      temp = nose_lefteye[i][j] + nose_righteye[i][j] + nose_leftmouth[i][j] + nose_rightmouth[i][j] + nose_chin[i][j] + my_potentials[NOSE][i][j];
      if (temp > max)
      {
        max = temp;
        answer.coords[NOSE].row = i;
        answer.coords[NOSE].col = j;
        assignments_D[NOSE].row=i;
        assignments_D[NOSE].col=j;  //assignments_D and answer.coords are both same thing, one is used for draw config
      }
    }
  } 


  //assignments to coords of GTImage that will be returned. As you can see they are acquired by using assignment of root (left eye) as index and
  //using it on appropriate coor_map
  answer.coords[LEFT_EYE].row = coor_map_n_le[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[LEFT_EYE].col = coor_map_n_le[answer.coords[NOSE].row][answer.coords[NOSE].col].col;
 
  answer.coords[RIGHT_EYE].row = coor_map_n_re[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[RIGHT_EYE].col = coor_map_n_re[answer.coords[NOSE].row][answer.coords[NOSE].col].col;

  answer.coords[LEFT_MOUTH].row = coor_map_n_lm[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[LEFT_MOUTH].col = coor_map_n_lm[answer.coords[NOSE].row][answer.coords[NOSE].col].col;

  answer.coords[RIGHT_MOUTH].row = coor_map_n_rm[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[RIGHT_MOUTH].col = coor_map_n_rm[answer.coords[NOSE].row][answer.coords[NOSE].col].col;

  answer.coords[CHIN].row = coor_map_n_c[answer.coords[NOSE].row][answer.coords[NOSE].col].row;
  answer.coords[CHIN].col = coor_map_n_c[answer.coords[NOSE].row][answer.coords[NOSE].col].col;

  assignments_D[LEFT_EYE].row = coor_map_n_le[assignments_D[NOSE].row][assignments_D[NOSE].col].row;
  assignments_D[LEFT_EYE].col = coor_map_n_le[assignments_D[NOSE].row][assignments_D[NOSE].col].col;
 
  assignments_D[RIGHT_EYE].row = coor_map_n_re[assignments_D[NOSE].row][assignments_D[NOSE].col].row;
  assignments_D[RIGHT_EYE].col = coor_map_n_re[assignments_D[NOSE].row][assignments_D[NOSE].col].col;

  assignments_D[LEFT_MOUTH].row = coor_map_n_lm[assignments_D[NOSE].row][assignments_D[NOSE].col].row;
  assignments_D[LEFT_MOUTH].col = coor_map_n_lm[assignments_D[NOSE].row][assignments_D[NOSE].col].col;

  assignments_D[RIGHT_MOUTH].row = coor_map_n_rm[assignments_D[NOSE].row][assignments_D[NOSE].col].row;
  assignments_D[RIGHT_MOUTH].col = coor_map_n_rm[assignments_D[NOSE].row][assignments_D[NOSE].col].col;

  assignments_D[CHIN].row = coor_map_n_c[assignments_D[NOSE].row][assignments_D[NOSE].col].row;
  assignments_D[CHIN].col = coor_map_n_c[assignments_D[NOSE].row][assignments_D[NOSE].col].col;
  
  
  int strtIndx=filename.find("_"), endIndx=filename.find(".");
  string imageNumber=filename.substr(strtIndx+1,endIndx-strtIndx-1);
  string imgOutput="D_Output_" + imageNumber + ".png";
  draw_configuration(imgOutput.c_str(), filename.c_str(), assignments_D);
  return answer;
   
} 



/************************************************************************************************
                                     MODEL E (Loopy Belief) 
************************************************************************************************ */
/*
This method is used to set a default value for a matrix of type SDoublePlane
*/

void set_temp (SDoublePlane &t,int rows,int cols)
{
  for(int i = 0; i<rows ; i++)
  {
    for(int j =0; j<cols ; j++)
    {
      t[i][j] = 0;
    }
  }
}

/*This function calculates the phi(Li,Lj) for the loopy function - similar to the one used before
except that there is no corresponding coordinate assignment done here.
*/

SDoublePlane Phi_calculator_loopy(Coordinate mean, Coordinate var,SDoublePlane my_potential) //returns an array of coordinates, which is the best possible assignments for two parts i and j
{
  int rows = my_potential.rows();
  int cols = my_potential.cols();
  SDoublePlane phi = SDoublePlane(rows,cols);
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      phi[i][j] = -1e100;
    }

  }
  for(int ix = 0; ix<rows; ix++)
  {
    for(int iy = 0; iy<cols; iy++)
    {
      double max_phi = -1e100;
      for (int jx = 0; jx<rows; jx+=STEP_SIZE)
      {
        for (int jy = 0; jy<cols; jy+=STEP_SIZE)
        {
          double t1 = pow(ix-jx-mean.row,2)/(2*var.row);
          double t2 = pow(iy-jy-mean.col,2)/(2*var.col);
          double temp = -t1-t2+my_potential[jx][jy];
          if(temp>max_phi)
          {
            max_phi = temp;
          }
        }
      }
      phi[ix][iy] = max_phi;
    }
  }
  return phi;

} 

/*
Function that takes in 2 SDoublePlane matrices a and b and returns another SDoublePlane that is the sum
of the corresponding elements in a and b.
*/

SDoublePlane add (SDoublePlane a,SDoublePlane b,int rows,int cols)
{
  SDoublePlane result = SDoublePlane(rows,cols);
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      result[i][j] = -1e100;
    }

  }
  for (int i = 0; i < rows ; i++ )
  {
    for ( int j = 0; j < cols ; j++ )
    {
      result[i][j] = a[i][j] + b[i][j];
    }
  }
  return result;
}

/*
Function which takes two SDoublePlane matrices and assigns the rhs to the lhs.
*/

void assign (SDoublePlane &lhs,SDoublePlane rhs,int rows,int cols)
{
  for (int i =0; i<rows;i++)
  {
    for (int j=0;j<cols;j++)
    {
      lhs[i][j] = rhs[i][j];
    }
  }
}

/*
Main Loopy Belief propagation method for model E of the assignment
*/

GTImage loopy(vector <SDoublePlane>my_potentials, string filename, Coordinate *means[6], Coordinate *vars[6])
{
  /*
  Here each node of the graph maintains two tables, one which stores the updated belief state at time
  t+1 and another for its belief state at time t.
  We denote them as node_message (for time t+1) and node_message_old (for time t).
  */
  //This part is used to make sure that means for part pairs (Li,Lj) and (Lj,Li) are different 
  // and negative of each other.
  for(int i=0; i<PART_COUNT; i++) 
  {
    for(int j=0; j<i; j++)
    {
        means[i][j].row = - means[j][i].row;
        means[i][j].col = - means[j][i].col;
        vars[i][j].row = vars[j][i].row;
        vars[i][j].col = vars[j][i].col;
    }
  }

  cout<<"\n Starting Loopy..."<<endl;
  //getting details about the image used
  SDoublePlane image1 = SImageIO::read_png_file(filename.c_str());
  int rows = image1.rows();//number of rows in the image
  int cols = image1.cols(); //number of cols in the image
  

  
  //Initialize the messages
  //node_message is the updated node_table
  vector <SDoublePlane>node_message (6,SDoublePlane(rows,cols)); 
  //node_message_old is the node table at time t-1
  vector <SDoublePlane>node_message_old (6,SDoublePlane(rows,cols));
  //give a default value to all the messages
  for ( int ii=0 ; ii<6 ; ii++ )
  {
    for ( int i=0 ; i<rows ; i++ )
    {
      for ( int j=0 ; j<cols ; j++ )
      {
        
          node_message[ii][i][j] = 0;
          node_message_old[ii][i][j] = 0;
      }
    }
  }
  //The loop starts here - i denotes the number of iterations of the loop
 

  for (int i = 0; i<1 ; i++)
  {
    cout<<"Iteration "<<i+1<<endl;
    //This loop runs for all the nodes in the graph. Each node updates the belief state of its neighbours
    //using its own belief state, the belief states of all its other neighnours and the phi between it and 
    //its neighbour.
    for (int node = 0; node < 6; node++ )
    {
      cout<<"Computing updates for node "<<node<<endl;
      if ( node == 0 )
      {
        // Update belief table for Left Mouth and Right Eye
        assign(node_message[3],add(add(node_message_old[3],add(add(node_message_old[1],my_potentials[0],rows,cols),Phi_calculator_loopy(means[0][3],vars[0][3],my_potentials[3]),rows,cols),rows,cols),node_message_old[0],rows,cols),rows,cols);
        assign(node_message[1],add(add(node_message_old[1],add(add(node_message_old[3],my_potentials[0],rows,cols),Phi_calculator_loopy(means[0][1],vars[0][1],my_potentials[1]),rows,cols),rows,cols),node_message_old[0],rows,cols),rows,cols);
        
      }
      if ( node == 1 )
      {
        //Update belief table for Left Eye, Nose and Right Mouth
        assign(node_message[0] , add(add(add(node_message_old[0],node_message_old[2],rows,cols),add(add(node_message_old[4],my_potentials[1],rows,cols),Phi_calculator_loopy(means[1][0],vars[1][0],my_potentials[0]),rows,cols),rows,cols),node_message_old[1],rows,cols),rows,cols);
        assign(node_message[2] , add(add(add(node_message_old[2],node_message_old[0],rows,cols),add(add(node_message_old[4],my_potentials[1],rows,cols),Phi_calculator_loopy(means[1][2],vars[1][2],my_potentials[2]),rows,cols),rows,cols),node_message_old[1],rows,cols),rows,cols);
        assign(node_message[4] , add(add(add(node_message_old[4],node_message_old[0],rows,cols),add(add(node_message_old[2],my_potentials[1],rows,cols),Phi_calculator_loopy(means[1][4],vars[1][4],my_potentials[4]),rows,cols),rows,cols),node_message_old[1],rows,cols),rows,cols);
      }
      if ( node == 2 )
      {
        //Update belief table for Right Eye and Chin
        assign(node_message[1] , add(add(add(node_message_old[1],node_message_old[5],rows,cols),add(my_potentials[2],Phi_calculator_loopy(means[2][1],vars[2][1],my_potentials[1]),rows,cols),rows,cols),node_message_old[2],rows,cols),rows,cols);
        assign(node_message[5] , add(add(add(node_message_old[5],node_message_old[1],rows,cols),add(my_potentials[2],Phi_calculator_loopy(means[2][5],vars[2][5],my_potentials[5]),rows,cols),rows,cols),node_message_old[2],rows,cols),rows,cols);
      }
      if ( node == 3 )
      {
        //Update belief table for Left Eye and Chin
        assign(node_message[0] , add(add(add(node_message_old[0],node_message_old[5],rows,cols),add(my_potentials[3],Phi_calculator_loopy(means[3][0],vars[3][0],my_potentials[0]),rows,cols),rows,cols),node_message_old[3],rows,cols),rows,cols);
        assign(node_message[5] , add(add(add(node_message_old[5],node_message_old[0],rows,cols),add(my_potentials[3],Phi_calculator_loopy(means[3][5],vars[3][5],my_potentials[5]),rows,cols),rows,cols),node_message_old[3],rows,cols),rows,cols);
      }
      if ( node == 4)
      {
        //Update belief table for Right Eye and Chin
        assign(node_message[1] , add(add(add(node_message_old[1],node_message_old[5],rows,cols),add(my_potentials[4],Phi_calculator_loopy(means[4][1],vars[4][1],my_potentials[1]),rows,cols),rows,cols),node_message_old[4],rows,cols),rows,cols);
        assign(node_message[5] , add(add(add(node_message_old[5],node_message_old[1],rows,cols),add(my_potentials[4],Phi_calculator_loopy(means[4][5],vars[4][5],my_potentials[5]),rows,cols),rows,cols),node_message_old[4],rows,cols),rows,cols);
      }
      if( node == 5 )
      {
        //Update belief table for Nose, Left Mouth and Right Mouth
        assign(node_message[2] , add(add(add(node_message_old[2],node_message_old[3],rows,cols),add(add(node_message_old[4],my_potentials[5],rows,cols),Phi_calculator_loopy(means[5][2],vars[5][2],my_potentials[2]),rows,cols),rows,cols),node_message_old[5],rows,cols),rows,cols);
        assign(node_message[3] , add(add(add(node_message_old[3],node_message_old[2],rows,cols),add(add(node_message_old[4],my_potentials[5],rows,cols),Phi_calculator_loopy(means[5][3],vars[5][3],my_potentials[3]),rows,cols),rows,cols),node_message_old[5],rows,cols),rows,cols);
        assign(node_message[4] , add(add(add(node_message_old[4],node_message_old[2],rows,cols),add(add(node_message_old[3],my_potentials[5],rows,cols),Phi_calculator_loopy(means[5][4],vars[5][4],my_potentials[4]),rows,cols),rows,cols),node_message_old[5],rows,cols),rows,cols);
      }
    }//end-for for all nodes
    //Assign the new node tables as the old tables.
    for( int n = 0; n<6 ; n++ )
    {
      assign(node_message_old[n],node_message[n],rows,cols);
    }
  }//end-for main loop
  //This part finds the best assignment of pixels for each node using its updated belief state.
  vector <Coordinate>best_configuration(6);
  for ( int node = 0; node < 6 ; node++ )
  {
    double temp = -1e100; //keeps a track of the maximum pixel value
    Coordinate temp_cord;

    for (int i = 0; i<rows ; i++)
    {
      for ( int j = 0; j<cols ; j++)
      {
        if (node_message[node][i][j] >= temp) //Check if the value is greater than the previous max value
        {
          temp = node_message[node][i][j];
          temp_cord.row = i;
          temp_cord.col = j;
        }
      }
    }
    best_configuration[node].row = temp_cord.row;
    cout<<"Node: "<<node<<" Row: "<<best_configuration[node].row<<endl;
    best_configuration[node].col = temp_cord.col;
    cout<<"Node: "<<node<<" Col: "<<best_configuration[node].col<<endl;
  }

  int strtIndx=filename.find("_"), endIndx=filename.find(".");
  string imageNumber=filename.substr(strtIndx+1,endIndx-strtIndx-1);
  string imgOutput="Loopy_Output_" + imageNumber + ".png";
  draw_configuration(imgOutput.c_str(), filename.c_str(), best_configuration);
  //draw_configuration("LOOPY.png", filename.c_str(), best_configuration);
  
  GTImage answer;
  answer.image_filename = filename;
  for (int i=0; i<6; i++)
  {
    answer.coords[i].row = best_configuration[i].row;  
    answer.coords[i].col = best_configuration[i].col;
  }
  return answer;
  
  

} //end of method

/*
This is heart of our program. This calls all the 4 inferences techniquies on each image in the testing data set. 
This also computes the accuracy of every inference method.
This is computed as:
 Total number of correctly placed parts / Total number of images
 where a part is classified as correctly placed if:
  | X1 - X2 | + | Y1 +Y2 | <= 20
  Where:
    X1 = X coordinate of part as computed
    X2 = X coordinate of part as per ground truth
    Y1 = Y coordinate of part as computed
    Y2 = Y coordinate of part as per ground truth
*/

void evaluation(vector<GTImage> testing_data_truth, vector<GTImage> training_data)
{
  
  
  UnaryPotential unary_potentials(training_data, 15);

  Coordinate *means[6]; //parameters for training
  Coordinate *vars[6];
  for (int i=0; i<6; i++)
  {
    means[i] = new Coordinate[6];
    vars[i] = new Coordinate[6];
  }

  //call learning
  learning(training_data, means, vars);
  //declare 4 accuracy variables and counter variables for 4 graph models B,C,D,E
  double B_accuracy;
  double C_accuracy;
  double D_accuracy;
  double E_accuracy;
  double B_counter=0;
  double C_counter=0;
  double D_counter=0;
  double E_counter=0;
  double total_parts = testing_data_truth.size() * PART_COUNT;

  //Iterating through all the images in the testing data set.
  for (int i=0; i< testing_data_truth.size(); i++)   
  {
    vector<SDoublePlane> my_potentials = unary_potentials.compute_unary_logpotential(testing_data_truth[i].image_filename);
    cout << "Running all 4 Inferences for file " << testing_data_truth[i].image_filename << endl;
    //Perform all inference methods for the current file
    GTImage temp_B = MAP_Inference_B(my_potentials, testing_data_truth[i].image_filename, means, vars);
    GTImage temp_C = MAP_Inference_C(my_potentials, testing_data_truth[i].image_filename, means, vars);
    GTImage temp_D = MAP_Inference_D(my_potentials, testing_data_truth[i].image_filename, means, vars);
    GTImage temp_E = loopy(my_potentials, testing_data_truth[i].image_filename, means, vars);

    // Compute number of correctly placed parts in each model
    for (int k=0; k<6; k++)
    {
      
      
      if ( ( abs(temp_B.coords[k].row - testing_data_truth[i].coords[k].row) <= 20) && (abs( temp_B.coords[k].col - testing_data_truth[i].coords[k].col) ) <= 20 )
      {
        B_counter++;
        
      }

      if ( ( abs(temp_C.coords[k].row - testing_data_truth[i].coords[k].row) <= 20) && (abs( temp_C.coords[k].col - testing_data_truth[i].coords[k].col) ) <= 20 )
      {
        
        C_counter++;
        
      } 

      if ( ( abs(temp_D.coords[k].row - testing_data_truth[i].coords[k].row) <= 20) && (abs( temp_D.coords[k].col - testing_data_truth[i].coords[k].col) ) <= 20 )
      {
        D_counter++;
        
      } 

      if ( ( abs(temp_E.coords[k].row - testing_data_truth[i].coords[k].row) <=20) && (abs( temp_E.coords[k].col - testing_data_truth[i].coords[k].col) ) <= 20 )
      {
        E_counter++; 
      } 

    } 
  }

  B_accuracy = B_counter /  total_parts;
  C_accuracy = C_counter /  total_parts;
  D_accuracy = D_counter /  total_parts;
  E_accuracy = E_counter / total_parts; 
  

  cout << "\n********************* RESULTS **************************" << endl;
  cout << "Accuracy of MAP Max Product for Model B = " << B_accuracy << endl;
  cout << "Accuracy of MAP Max Product for Model C = " << C_accuracy << endl;
  cout << "Accuracy of MAP Max Product for Model D = " << D_accuracy << endl;
  cout << "Accuracy of MAP Max Product for Model E = " << E_accuracy << endl; 


} 






int main(int argc, char *argv[])
{
  try {
    if(!(argc == 3)) {
      cerr << "usage: " << argv[0] << " training_gt_file testing_gt_file" << endl;
      return 1;
    }

    string train_gt_filename(argv[1]), test_gt_filename(argv[2]);

    // This is a parameter used to compute the unary potentials. It specifies
    //  how big the part appearance models should be.
    const int template_size=15;

    // Load in the training data
    vector<GTImage> training_data = load_groundtruth(train_gt_filename);
    vector<GTImage> testing_data = load_groundtruth(test_gt_filename);

    // Now set up the unary potential functions
    UnaryPotential unary_potentials(training_data, template_size);
    cerr << "Done learning unary potentials" << endl << endl;

    // Now run naive inference on one of the sample images
    //  
    string sample_filename = "imgs/image_0091.png";
    cerr << "Running naive inference on " << sample_filename << "..." << endl;

    // Compute unary potentials for this image
    vector<SDoublePlane> my_potentials = unary_potentials.compute_unary_logpotential(sample_filename);
    // Naive inference just involves choosing best coordinate in each individual
    //  unary potential function.
    const int part_count=6;
    vector<Coordinate> best_configuration(part_count);
    for(int ii=0; ii<part_count; ii++)
      {
  double max_val=-1e100;
  for(int i=0; i<my_potentials[ii].rows(); i++)
    for(int j=0; j<my_potentials[ii].cols(); j++)
      if(my_potentials[ii][i][j] > max_val)
        {
    best_configuration[ii].row = i;
    best_configuration[ii].col = j;
    max_val = my_potentials[ii][i][j];
    //cout << my_potentials[ii][i][j] <<endl;
    //cout << "max" << max_val << endl;
        }
      }

    // Visualize result by drawing detected parts on the image, and 
    //  outputing to a new image file.
    //
    draw_configuration("naive_result.png", sample_filename.c_str(), best_configuration);
    //Calling evaluation function which performs inference for all images in the testing data set for all models from b-e.



  Coordinate *means[6]; //parameters for training
  Coordinate *vars[6];
  for (int i=0; i<6; i++)
  {
    means[i] = new Coordinate[6];
    vars[i] = new Coordinate[6];
  }

  //call learning
  learning(training_data, means, vars);
  evaluation(testing_data, training_data);

  //MAP_Inference_B(my_potentials, sample_filename, means, vars);
  //MAP_Inference_C(my_potentials, sample_filename, means, vars);
  //MAP_Inference_D(my_potentials, sample_filename, means, vars);

    

  } catch(std::string &str) {
    cerr << str << endl;
  }

  return 0;
} 
