#pragma once

#include "Utils.h"
#include "Matrix.h"
#include "Defs.h"


void ToHessenbergForm(Matrix &m);

void IterationQR(Matrix &h);

std::vector<Complex> EigQR(Matrix m, int &iter = it);

