//
//  volatility.hpp
//  pde_solver
//
//  Created by John Sibony on 25/01/2019.
//  Copyright © 2019 John Sibony. All rights reserved.
//

#ifndef volatility_hpp
#define volatility_hpp


class Volatility
{
    
public:
    
    Volatility();
    Volatility(double vol);
    double vol() const;
    
private:
    
    const double _vol;

};

#endif /* volatility_hpp */

