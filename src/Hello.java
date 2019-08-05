/*
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/
/*
 * Author of revised version: Franklyn Pinedo
 *
 * Adapted from Source Code in C of Tutorial/User's Guide for MPI by
 * Peter Pacheco.
 */
/*
 * Copyright (c) 2011      Cisco Systems, Inc.  All rights reserved.
 *
 */

import mpi.*;
import java.nio.*;

class Hello {

    params inputParams = new params();

    static public void main(String[] args) throws MPIException {


		MPI.Init(args);

		int myrank = MPI.COMM_WORLD.getRank();
		int size = MPI.COMM_WORLD.getSize() ;

		params inputParams = new params();



		System.out.println("Hello world from rank "  + " of " + size);
        inputParams.print();

		MPI.Finalize();
	}
}
