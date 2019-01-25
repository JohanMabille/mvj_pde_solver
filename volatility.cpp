//
//  volatility.cpp
//  pde_solver
//
//  Created by John Sibony on 25/01/2019.
//  Copyright © 2019 John Sibony. All rights reserved.
//

#include "volatility.hpp"

Volatility::Volatility(): _vol(0.2)
{
    
}

Volatility::Volatility(double vol): _vol(vol)
{
    
}

double Volatility::vol() const
{
    return _vol;
}
