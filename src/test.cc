#include <vector>
#include <iostream>
#include <random>

#include "../include/hmm.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Usage: ./baum_welch <initialization file>" << endl;
        return 0;
    }    

    HMM* model = new HMM();
    model->initialize_vectors(argv[1]);

    cout << "pi: " << endl;
    for (int i=0; i<model->M; i++)
        cout << model->pi[i] << " ";
    
    cout << endl;

    cout << "A: " << endl;
    for (int i=0; i<model->M; i++) {
        for (int j=0; j<model->M; j++)
            cout << model->A[i][j] << " ";
        cout << endl;
    }

    cout << "B: " << endl;
    for (int i=0; i<model->M; i++) {
        for (int j=0; j<model->N; j++)
            cout << model->B[i][j] << " ";
        cout << endl;
    }
    
    return 1;

}
