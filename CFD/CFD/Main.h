#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Windows.h>
#include <iostream>
#include <string.h>
#include "elliptic.h"
#include "hyperbolic.h"
#include "parabolic.h"
#include "TDMA.h"
int Eqnselector(void);
void starter(int Eqn);
void repeater(void);