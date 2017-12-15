#pragma once
/*////////////////////////////////////////////////////////////////////////////
WingGenerators.h

Methods to generate wing shapes.

Copyright 2017 HJA Bird

mFlow is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mFlow is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mFlow.  If not, see <http://www.gnu.org/licenses/>.
*/////////////////////////////////////////////////////////////////////////////

#include "WingProjectionGeometry.h"

namespace mFlow {
	namespace WingGenerators {

		// Set wing_obj to an elliptically shaped wing of given span and aspect_ratio
		void elliptic(WingProjectionGeometry &wing_obj, double span, double aspect_ratio);
		// Set wing_obj to an rectangularly shaped wing of given span and aspect_ratio
		void rectangular(WingProjectionGeometry &wing_obj, double span, double aspect_ratio);

	}
}
