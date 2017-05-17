#define _USE_MATH_DEFINES
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cmath>

#include "classifier.h"
#include "EasyBMP.h"
#include "linear.h"
#include "argvparser.h"
#include "matrix.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;

using CommandLineProcessing::ArgvParser;

typedef vector<pair<BMP*, int> > TDataSet;
typedef vector<pair<string, int> > TFileList;
typedef vector<pair<vector<float>, int> > TFeatures;


// Load list of files and its labels from 'data_file' and
// stores it in 'file_list'
void LoadFileList(const string& data_file, TFileList* file_list) {
    ifstream stream(data_file.c_str());

    string filename;
    int label;
    
    int char_idx = data_file.size() - 1;
    for (; char_idx >= 0; --char_idx)
        if (data_file[char_idx] == '/' || data_file[char_idx] == '\\')
            break;
    string data_path = data_file.substr(0,char_idx+1);
    
    while(!stream.eof() && !stream.fail()) {
        stream >> filename >> label;
        if (filename.size())
            file_list->push_back(make_pair(data_path + filename, label));
    }

    stream.close();
}

// Load images by list of files 'file_list' and store them in 'data_set'
void LoadImages(const TFileList& file_list, TDataSet* data_set) {
    for (size_t img_idx = 0; img_idx < file_list.size(); ++img_idx) {
            // Create image
        BMP* image = new BMP();
            // Read image from file
        image->ReadFromFile(file_list[img_idx].first.c_str());
            // Add image and it's label to dataset
        data_set->push_back(make_pair(image, file_list[img_idx].second));
    }
}

// Save result of prediction to file
void SavePredictions(const TFileList& file_list,
                     const TLabels& labels, 
                     const string& prediction_file) {
        // Check that list of files and list of labels has equal size 
    assert(file_list.size() == labels.size());
        // Open 'prediction_file' for writing
    ofstream stream(prediction_file.c_str());

        // Write file names and labels to stream
    for (size_t image_idx = 0; image_idx < file_list.size(); ++image_idx)
        stream << file_list[image_idx].first << " " << labels[image_idx] << endl;
    stream.close();
}

Matrix<int> makeGrey(BMP* img){ //creating grey image
    Matrix<int> greyImg = Matrix<int>(img->TellWidth(), img->TellHeight());
        for (int i = 0; i < img->TellWidth(); ++i) {
            for (int j = 0; j < img->TellHeight(); ++j)
            {
                RGBApixel pixel = img->GetPixel(i, j);
                greyImg(i, j) = 0.299 * pixel.Red + 0.587 * pixel.Green + 0.114 * pixel.Blue;
            }
        }
    return greyImg;
}

Matrix<float> makeTans(Matrix<int>& in1, Matrix<int>& in2){
    Matrix<float> tanMatrix =  Matrix<float>(in1.n_rows, in1.n_cols);
    for (uint i = 0; i < in1.n_rows; ++i) {
        for (uint j = 0; j < in1.n_cols; ++j) {    
            int pixel1 = in1(i,j);
            int pixel2 = in2(i,j);

            tanMatrix(i,j) = atan2(pixel2,pixel1);
        }
    }
    return tanMatrix;
}

Matrix<int> makeMods(Matrix<int> &in1, Matrix<int> &in2){
    Matrix<int> modMatrix = Matrix<int>(in1.n_rows, in1.n_cols);
    for (uint i = 0; i < in1.n_rows; ++i) {
        for (uint j = 0; j < in1.n_cols; ++j) {    
            int pixel1 = in1(i,j);
            int pixel2 = in2(i,j);

            modMatrix(i,j) = sqrt(pixel2 * pixel2 + pixel1 * pixel1);
        }
    }
    return modMatrix;
}

vector<float> hystogramCreator(Matrix<int> mods, Matrix<float> tans, unsigned int segTotal) {
    vector<float> hystogram(segTotal, 0.0);
    float segSize = 2 * M_PI / segTotal;
    for (uint i = 0; i < mods.n_rows; ++i) {
        for (uint j = 0; j < mods.n_cols; ++j) {
            uint hystIdx = (tans(i, j) + M_PI) / segSize;
            hystogram[hystIdx] += mods(i, j);
        }
    }
    float sqrs = 0;
    for (uint i = 0; i < segTotal; ++i){
        sqrs += hystogram[i] * hystogram[i];
    }
    float eps = 0.003;
    float norm = sqrt(sqrt(sqrs) + eps*eps);
    for (uint i = 0; i < segTotal; ++i){
        hystogram[i] /= norm;
    }
    return hystogram;
}

vector<float> hystogramLbp(Matrix<int> in, uint segTotal) {
    vector<float> hystogram(segTotal, 0.0);
    for (uint i = 0; i < in.n_rows; ++i) {
        for (uint j = 0; j < in.n_cols; ++j) {
            hystogram[in(i, j)] += 1;
        }
    }
    float sqrs = 0;
    for (uint i = 0; i < segTotal; ++i){
        sqrs += hystogram[i] * hystogram[i];
    }
    float eps = 0.003;
    float norm = sqrt(sqrt(sqrs) + eps*eps);
    for (uint i = 0; i < segTotal; ++i){
        hystogram[i] /= norm;
    }
    return hystogram;
}

Matrix<RGBApixel> setAvgColors(Matrix<RGBApixel> img){
    Matrix<RGBApixel> avgColors(img.n_rows, img.n_cols);
    uint sumRed = 0, sumGreen = 0, sumBlue = 0;
    for (uint i = 0; i < img.n_rows; ++i) {
        for (uint j = 0; j < img.n_cols; ++j)
        {
            sumRed += img(i,j).Red;
            sumGreen += img(i,j).Green;
            sumBlue += img(i,j).Blue;
        }
    }
    for (uint i = 0; i < avgColors.n_rows; ++i) {
        for (uint j = 0; j < avgColors.n_cols; ++j)
        {
            avgColors(i,j).Red = sumRed/(avgColors.n_rows * avgColors.n_cols);
            avgColors(i,j).Green = sumGreen/(avgColors.n_rows * avgColors.n_cols);
            avgColors(i,j).Blue = sumBlue/(avgColors.n_rows * avgColors.n_cols);
        }
    }
    return avgColors;
}

// Extract features from dataset.
// You should implement this function by yourself =)
void ExtractFeatures(const TDataSet& data_set, TFeatures* features) {
    for (size_t image_idx = 0; image_idx < data_set.size(); ++image_idx) {
        
        // PLACE YOUR CODE HERE
        BMP* img = data_set[image_idx].first; // img extraction
        Matrix<int> sobelHorKern = {
                        {0, 0, 0},
                        {-1, 0, 1},
                        {0, 0, 0} };
        Matrix<int> sobelVertKern = {
                        {0, 1, 0},
                        {0, 0, 0},
                        {0, -1, 0} };
        Matrix<int> lbpKernel = {
                        {128, 64, 32},
                        {1, 0, 16},
                        {2, 4, 8} };
        Filter<int> sobelHor(sobelHorKern, 1, 1), sobelVert(sobelVertKern, 1, 1);
        Filter<int> lbpUser(lbpKernel, 1, 1, 1);
        //creating matrixes of atans & mods
        Matrix<int> greyImage = makeGrey(img);
        Matrix<int> horSobel = greyImage.unary_map(sobelHor);
        Matrix<int> vertSobel = greyImage.unary_map(sobelVert);
        Matrix<int> mods = makeMods(horSobel, vertSobel);
        Matrix<float> tans = makeTans(horSobel, vertSobel);
        
        //LBP matrix
        Matrix<int> lbpMatrix = greyImage.unary_map(lbpUser);
        uint lbpSegments = 256;
        vector< vector<float> > lbpCells;
        //hystogram step

        vector< vector<float> > cellHystograms;
        unsigned int hystXlen = 8;
        unsigned int hystYlen = 8;
        unsigned int segTotal = 32;
        uint usualMatXsize = mods.n_rows / hystXlen;
        uint usualMatYsize = mods.n_cols / hystYlen;
        for (uint i = 0; i < usualMatXsize * hystXlen; i += usualMatXsize) {
            for (uint j = 0; j < usualMatYsize * hystYlen; j+= usualMatYsize) {
                vector<float> hystogram = hystogramCreator(mods.submatrix(i, j, usualMatXsize, usualMatYsize),
                tans.submatrix(i, j, usualMatXsize, usualMatYsize), segTotal);

                vector<float> lbpHystogram = hystogramLbp(lbpMatrix.submatrix(i, j, usualMatXsize, usualMatYsize), lbpSegments);
                lbpCells.push_back(lbpHystogram); // lbp Hystos
                cellHystograms.push_back(hystogram);
            }
        }
        vector<float> one_image_features;
        for (auto element = cellHystograms.begin(); element != cellHystograms.end(); ++element){
            one_image_features.insert(one_image_features.cbegin(), (*element).cbegin(), (*element).cend());
        }
        // lbp hystos append
        for (auto element = lbpCells.begin(); element != lbpCells.end(); ++element){
            one_image_features.insert(one_image_features.cbegin(), (*element).cbegin(), (*element).cend());
        }

        Matrix<RGBApixel> imgMatrix(img->TellWidth(), img->TellHeight());

        uint xPicBlockNum = 8;
        uint yPicBlockNum = 8;
        usualMatXsize = imgMatrix.n_rows / xPicBlockNum;
        usualMatYsize = imgMatrix.n_cols / yPicBlockNum;
        vector< Matrix<RGBApixel> > matVector;

        for (uint i = 0; i < usualMatXsize * xPicBlockNum; i += usualMatXsize) {
            for (uint j = 0; j < usualMatYsize * yPicBlockNum; j+= usualMatYsize) {
                Matrix<RGBApixel> avgColors = setAvgColors(imgMatrix.submatrix(i, j, usualMatXsize, usualMatYsize));
                matVector.push_back(avgColors);
            }
        }
        vector <float> normColorMap;

        for (auto element = matVector.begin(); element != matVector.end(); ++element){
            normColorMap.push_back( ((*element)(0,0).Red / 255.0) );
            normColorMap.push_back( ((*element)(0,0).Green / 255.0) );
            normColorMap.push_back( ((*element)(0,0).Blue / 255.0) );
        }

        one_image_features.insert(one_image_features.cbegin(), normColorMap.cbegin(), normColorMap.cend());
        
        features->push_back(make_pair(one_image_features, data_set[image_idx].second));
    }
}

// Clear dataset structure
void ClearDataset(TDataSet* data_set) {
        // Delete all images from dataset
    for (size_t image_idx = 0; image_idx < data_set->size(); ++image_idx)
        delete (*data_set)[image_idx].first;
        // Clear dataset
    data_set->clear();
}

// Train SVM classifier using data from 'data_file' and save trained model
// to 'model_file'
void TrainClassifier(const string& data_file, const string& model_file) {
        // List of image file names and its labels
    TFileList file_list;
    
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // Model which would be trained
    TModel model;
        // Parameters of classifier
    TClassifierParams params;
    
        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // PLACE YOUR CODE HERE
        // You can change parameters of classifier here
    params.C = 0.01;
    TClassifier classifier(params);
        // Train classifier
    classifier.Train(features, &model);
        // Save model to file
    model.Save(model_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

// Predict data from 'data_file' using model from 'model_file' and
// save predictions to 'prediction_file'
void PredictData(const string& data_file,
                 const string& model_file,
                 const string& prediction_file) {
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // List of image labels
    TLabels labels;

        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // Classifier 
    TClassifier classifier = TClassifier(TClassifierParams());
        // Trained model
    TModel model;
        // Load model from file
    model.Load(model_file);
        // Predict images by its features using 'model' and store predictions
        // to 'labels'
    classifier.Predict(features, model, &labels);

        // Save predictions
    SavePredictions(file_list, labels, prediction_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

int main(int argc, char** argv) {
    // Command line options parser
    ArgvParser cmd;
        // Description of program
    cmd.setIntroductoryDescription("Machine graphics course, task 2. CMC MSU, 2014.");
        // Add help option
    cmd.setHelpOption("h", "help", "Print this help message");
        // Add other options
    cmd.defineOption("data_set", "File with dataset",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("model", "Path to file to save or load model",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("predicted_labels", "Path to file to save prediction results",
        ArgvParser::OptionRequiresValue);
    cmd.defineOption("train", "Train classifier");
    cmd.defineOption("predict", "Predict dataset");
        
        // Add options aliases
    cmd.defineOptionAlternative("data_set", "d");
    cmd.defineOptionAlternative("model", "m");
    cmd.defineOptionAlternative("predicted_labels", "l");
    cmd.defineOptionAlternative("train", "t");
    cmd.defineOptionAlternative("predict", "p");

        // Parse options
    int result = cmd.parse(argc, argv);

        // Check for errors or help option
    if (result) {
        cout << cmd.parseErrorDescription(result) << endl;
        return result;
    }

        // Get values 
    string data_file = cmd.optionValue("data_set");
    string model_file = cmd.optionValue("model");
    bool train = cmd.foundOption("train");
    bool predict = cmd.foundOption("predict");

        // If we need to train classifier
    if (train)
        TrainClassifier(data_file, model_file);
        // If we need to predict data
    if (predict) {
            // You must declare file to save images
        if (!cmd.foundOption("predicted_labels")) {
            cerr << "Error! Option --predicted_labels not found!" << endl;
            return 1;
        }
            // File to save predictions
        string prediction_file = cmd.optionValue("predicted_labels");
            // Predict data
        PredictData(data_file, model_file, prediction_file);
    }
}