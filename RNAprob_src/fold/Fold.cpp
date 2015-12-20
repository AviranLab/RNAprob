/*
 * A program that folds a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "Fold.h"
#include <time.h>
///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
Fold::Fold() {

	// Initalize the calculation type description.
	calcType = "Single strand folding";

	// Initialize the "experimental" offset.
	experimentalOffset = 0.0;

	// Initialize the "experimental" scaling.
	experimentalScaling = 1.0;

	// Initialize the SHAPE intercept.
	intercept = -0.6;

	// Initialize the single stranded SHAPE intercept.
	interceptSingle = 0;

	// Initialize the nucleic acid type.
	isRNA = true;

	// Initialize the maximum internal bulge loop size.
	maxLoop = 30;

	// Initialize the maximum pairing distance between nucleotides.
	maxDistance = -1;

	// Initialize the maximum number of structures.
	maxStructures = 20;

	// Initialize the maximum percent energy difference.
	percent = 10;

	// Initialize the SHAPE slope.
	slope = 1.8;

	//Initialize the differntial SHAPE slope.
	Dslope = 2.11;

	// Initialize the single stranded SHAPE slope.
	slopeSingle = 0;

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the calculation temperature.
	bootstrap = 0;

	// Initialize the folding window size.
	windowSize = -1;

	//  Initialize the quickfold (mfe only) variable.
	quickfold = false;
	
	// FD
	// Initialize the bin size of the histogram
	binSize = 0.1;
	
	// Initialize the scheme to be three-state version
	twoStateVersion = false;
	
	// Initialize the decoder to be empirical version
	smoothVersion = false;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool Fold::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "Fold" );
	parser->addParameterDescription( "seq file", "The name of a file containing an input sequence." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Add the constraint file option.
	vector<string> constraintOptions;
	constraintOptions.push_back( "-c" );
	constraintOptions.push_back( "-C" );
	constraintOptions.push_back( "--constraint" );
	parser->addOptionFlagsWithParameters( constraintOptions, "Specify a constraints file to be applied. Default is to have no constraints applied." );

	// Add the DNA option.
	vector<string> dnaOptions;
	dnaOptions.push_back( "-d" );
	dnaOptions.push_back( "-D" );
	dnaOptions.push_back( "--DNA" );
	parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );
	
	// Add the twoState option.
	vector<string> twoStateOptions;
	twoStateOptions.push_back( "-2s" );
	twoStateOptions.push_back( "-2S" );
	parser->addOptionFlagsNoParameters( twoStateOptions, "Specify that scheme to be used is two-state version. Default is to use three-state version." );
	
	// Add the smooth version option.
	vector<string> smoothVersionOptions;
	smoothVersionOptions.push_back( "-smooth" );
	smoothVersionOptions.push_back( "-SMOOTH" );
	smoothVersionOptions.push_back( "--SMOOTH" );
	parser->addOptionFlagsNoParameters( smoothVersionOptions, "Specify that decoder to be smoothed version. Default is to use empirical version." );
	

	// Add the DNA option.
	vector<string> quickfoldOptions;
	quickfoldOptions.push_back( "-mfe" );
	quickfoldOptions.push_back( "-MFE" );
	quickfoldOptions.push_back( "--MFE" );
	parser->addOptionFlagsNoParameters( quickfoldOptions, "Specify that only the minimum free energy structure is needed.  No savefiles can be generated.  This saves nearly half the calculation time, but provides less information." );

	// Add the double stranded offset option.
	vector<string> doubleOffsetOptions;
	doubleOffsetOptions.push_back( "-dso" );
	doubleOffsetOptions.push_back( "-DSO" );
	doubleOffsetOptions.push_back( "--doubleOffset" );
	parser->addOptionFlagsWithParameters( doubleOffsetOptions, "Specify a double-stranded offset file, which adds energy bonuses to particular double-stranded nucleotides. Default is to have no double-stranded offset specified." );

	// Add the maximum loop size option.
	vector<string> loopOptions;
	loopOptions.push_back( "-l" );
	loopOptions.push_back( "-L" );
	loopOptions.push_back( "--loop" );
	parser->addOptionFlagsWithParameters( loopOptions, "Specify a maximum internal/bulge loop size. Default is 30 unpaired numcleotides." );

	// Add the maximum number of structures option.
	vector<string> maxStructuresOptions;
	maxStructuresOptions.push_back( "-m" );
	maxStructuresOptions.push_back( "-M" );
	maxStructuresOptions.push_back( "--maximum" );
	parser->addOptionFlagsWithParameters( maxStructuresOptions, "Specify a maximum number of structures. Default is 20 structures." );

	// Add the maximum pairing distance option.
	vector<string> distanceOptions;
	distanceOptions.push_back( "-md" );
	distanceOptions.push_back( "-MD" );
	distanceOptions.push_back( "--maxdistance" );
	parser->addOptionFlagsWithParameters( distanceOptions, "Specify a maximum pairing distance between nucleotides. Default is no restriction on distance between pairs." );

	// Add the percent energy difference option.
	vector<string> percentOptions;
	percentOptions.push_back( "-p" );
	percentOptions.push_back( "-P" );
	percentOptions.push_back( "--percent" );
	parser->addOptionFlagsWithParameters( percentOptions, "Specify a maximum percent energy difference. Default is 10 percent (specified as 10, not 0.1)." );

	// Add the save file option.
	vector<string> saveOptions;
	saveOptions.push_back( "-s" );
	saveOptions.push_back( "-S" );
	saveOptions.push_back( "--save" );
	parser->addOptionFlagsWithParameters( saveOptions, "Specify the name of a save file, needed for dotplots and refolding. Default is not to generate a save file." );

	// Add the SHAPE option.
	vector<string> shapeOptions;
	shapeOptions.push_back( "-sh" );
	shapeOptions.push_back( "-SH" );
	shapeOptions.push_back( "--SHAPE" );
	parser->addOptionFlagsWithParameters( shapeOptions, "Specify a SHAPE restraints file to be applied. These constraints are pseudoenergy restraints. Default is to have no restraints applied." );

	// Add the SHAPE option.
	vector<string> DshapeOptions;
	DshapeOptions.push_back( "-dsh" );
	DshapeOptions.push_back( "-DSH" );
	DshapeOptions.push_back( "--DSHAPE" );
	parser->addOptionFlagsWithParameters( DshapeOptions, "Specify a differential SHAPE restraints file to be applied. These constraints are pseudoenergy restraints. Default is to have no restraints applied." );

	// Add the DMS option.
	vector<string> dmsOptions;
	dmsOptions.push_back( "-dms" );
	dmsOptions.push_back( "-DMS" );
	dmsOptions.push_back( "--DMS" );
	parser->addOptionFlagsWithParameters( dmsOptions, "Specify a DMS restraints file to be applied. These restraints are pseudoenergy constraints. Default is to have no restraints applied." );

	// Add the CMCT option.
	vector<string> cmctOptions;
	cmctOptions.push_back( "-cmct" );
	cmctOptions.push_back( "-CMC" );
	cmctOptions.push_back( "--CMCT" );
	parser->addOptionFlagsWithParameters( cmctOptions, "Specify a CMCT constraints file to be applied. These constraints are pseudoenergy constraints. Default is to have no constraints applied." );


	// Add the SHAPE intercept option.
	vector<string> shapeInterceptOptions;
	shapeInterceptOptions.push_back( "-si" );
	shapeInterceptOptions.push_back( "-SI" );
	shapeInterceptOptions.push_back( "--SHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptOptions, "Specify an intercept used with SHAPE restraints. Default is -0.6 kcal/mol." );

	// Add the SHAPE slope option.
	vector<string> shapeSlopeOptions;
	shapeSlopeOptions.push_back( "-sm" );
	shapeSlopeOptions.push_back( "-SM" );
	shapeSlopeOptions.push_back( "--SHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeOptions, "Specify a slope used with SHAPE renstraints. Default is 1.8 kcal/mol." );

	// Add the Differentail SHAPE slope option.
	vector<string> DshapeSlopeOptions;
	DshapeSlopeOptions.push_back( "-dsm" );
	DshapeSlopeOptions.push_back( "-DSM" );
	DshapeSlopeOptions.push_back( "--DSHAPEslope" );
	parser->addOptionFlagsWithParameters( DshapeSlopeOptions, "Specify a slope used with differential SHAPE restraints. Default is 2.11 kcal/mol." );


	// Add the single stranded offset option.
	vector<string> singleOffsetOptions;
	singleOffsetOptions.push_back( "-sso" );
	singleOffsetOptions.push_back( "-SSO" );
	singleOffsetOptions.push_back( "--singleOffset" );
	parser->addOptionFlagsWithParameters( singleOffsetOptions, "Specify a single-stranded offset file, which adds energy bonuses to particular single-stranded nucleotides. Default is to have no single-stranded offset specified." );

	// Add the temperature option.
	vector<string> tempOptions;
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );

	// Add the bootstrap option.
	vector<string> bootstrapOptions;
	bootstrapOptions.push_back( "-boot" );
	bootstrapOptions.push_back( "-B" );
	bootstrapOptions.push_back( "--bootstrap" );
	parser->addOptionFlagsWithParameters( bootstrapOptions, "Specify the number of bootstrap iterations to be done to retrieve base pair confidence. Defaults to no bootstrapping." );


	// Add the unpaired SHAPE intercept option.
	vector<string> shapeInterceptUnpairedOptions;
	shapeInterceptUnpairedOptions.push_back( "-usi" );
	shapeInterceptUnpairedOptions.push_back( "-USI" );
	shapeInterceptUnpairedOptions.push_back( "--unpairedSHAPEintercept" );
	parser->addOptionFlagsWithParameters( shapeInterceptUnpairedOptions, "Specify an intercept used with unpaired SHAPE constraints. Default is 0 kcal/mol." );

	// Add the unpaired SHAPE slope option.
	vector<string> shapeSlopeUnpairedOptions;
	shapeSlopeUnpairedOptions.push_back( "-usm" );
	shapeSlopeUnpairedOptions.push_back( "-USM" );
	shapeSlopeUnpairedOptions.push_back( "--unpairedSHAPEslope" );
	parser->addOptionFlagsWithParameters( shapeSlopeUnpairedOptions, "Specify a slope used with unpaired SHAPE constraints. Default is 0 kcal/mol." );

	// Add the window size option.
	vector<string> windowOptions;
	windowOptions.push_back( "-w" );
	windowOptions.push_back( "-W" );
	windowOptions.push_back( "--window" );
	parser->addOptionFlagsWithParameters( windowOptions, "Specify a window size. Default is determined by the length of the sequence." );

	// Add the experimental pair bonus option.
	vector<string> experimentalOptions;
	experimentalOptions.push_back( "-x" );
	experimentalOptions.push_back( "-X" );
	experimentalOptions.push_back( "--experimentalPairBonus" );
	parser->addOptionFlagsWithParameters( experimentalOptions, "Input text file with bonuses (in kcal) as a matrix. As with SHAPE, bonuses will be applied twice to internal base pairs, once to edge base pairs, and not at all to single stranded regions. Default is no experimental pair bonus file specified." );

	// Add the experimental pair bonus offset option.
	vector<string> experimentalOffsetOptions;
	experimentalOffsetOptions.push_back( "-xo" );
	parser->addOptionFlagsWithParameters( experimentalOffsetOptions, "Specify an intercept (overall offset) to use with the 2D experimental pair bonus file. Default is 0.0 (no change to input bonuses)." );

	// Add the experimental pair bonus scaling option.
	vector<string> experimentalScalingOptions;
	experimentalScalingOptions.push_back( "-xs" );
	parser->addOptionFlagsWithParameters( experimentalScalingOptions, "Specify a number to multiply the experimental pair bonus matrix by. Default is 1.0 (no change to input bonuses)." );

	//FD
	//Add the bin size for the histogram
	vector<string> binSizeOptions;
	binSizeOptions.push_back( "-bs" );
	parser->addOptionFlagsWithParameters( binSizeOptions, "Specify the bin size for the histogram. Default is 0.1. " );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		seqFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Get the constraint file option.
	if( !parser->isError() ) { constraintFile = parser->getOptionString( constraintOptions, true ); }

	// Get the DNA option.
	if( !parser->isError() ) { isRNA = !parser->contains( dnaOptions ); }
	
	// FD
	// Get the twoState option.
	if( !parser->isError() ) { twoStateVersion = parser->contains( twoStateOptions ); }

	// FD
	// Get the smoothVersion option.
	if( !parser->isError() ) { smoothVersion = parser->contains( smoothVersionOptions ); }
	
	
	// Get the qickfold option.
	if( !parser->isError() ) { quickfold = parser->contains( quickfoldOptions ); }

	// Get the double stranded offset option.
	if( !parser->isError() ) { doubleOffsetFile = parser->getOptionString( doubleOffsetOptions, true ); }

	// Get the maximum loop size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( loopOptions, maxLoop );
		if( maxLoop < 0 ) { parser->setError( "maximum loop size" ); }
	}

	// Get the maximum distance option.
	if( !parser->isError() ) {
		parser->setOptionInteger( distanceOptions, maxDistance );
		bool badDistance =
		  ( maxDistance < 0 ) &&
		  ( maxDistance != -1 );
		if( badDistance ) { parser->setError( "maximum pairing distance" ); }
	}

	// Get the maximum number of structures option.
	if( !parser->isError() ) {
		parser->setOptionInteger( maxStructuresOptions, maxStructures );
		if( maxStructures <= 0 ) { parser->setError( "maximum number of structures" ); }
	}

	// Get the percent energy difference option.
	if( !parser->isError() ) {
		parser->setOptionDouble( percentOptions, percent );
		if( percent < 0 ) { parser->setError( "percent energy difference" ); }
	}

	// Get the save file option.
	if( !parser->isError() ) { saveFile = parser->getOptionString( saveOptions, false ); }

	// Set modifier type
	if( !parser->isError() ) {
		if(parser->contains(dmsOptions))
		    modifier = "DMS"; 
		if(parser->contains(shapeOptions))
		    modifier = "SHAPE"; 
		if(parser->contains(cmctOptions))
		    modifier = "CMCT"; 
	}

	// Get the SHAPE, DMS, and CMCT data and options.
	if( !parser->isError() ) {
		SHAPEFile = parser->getOptionString( shapeOptions );
		DSHAPEFile = parser->getOptionString( DshapeOptions );
		if( !parser->isError() ) { DMSFile = parser->getOptionString( dmsOptions ); }
		if( !parser->isError() ) { CMCTFile = parser->getOptionString( cmctOptions ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptOptions, intercept ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeOptions, slope ); }
		if( !parser->isError() ) { parser->setOptionDouble( DshapeSlopeOptions, Dslope ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeInterceptUnpairedOptions, interceptSingle ); }
		if( !parser->isError() ) { parser->setOptionDouble( shapeSlopeUnpairedOptions, slopeSingle ); }
	}

	// Get the single stranded offset option.
	if( !parser->isError() ) { singleOffsetFile = parser->getOptionString( singleOffsetOptions, true ); }

	// Get the temperature option.
	if( !parser->isError() ) {
		parser->setOptionDouble( tempOptions, temperature );
		if( temperature < 0 ) { parser->setError( "temperature" ); }
	}

	// Get the number of bootstraps.
	if( !parser->isError() ) {
		parser->setOptionDouble( bootstrapOptions, bootstrap );
		if( bootstrap < 0 ) { parser->setError( "bootstrap" ); }
	}

	// Get the window size option.
	if( !parser->isError() ) {
		parser->setOptionInteger( windowOptions, windowSize );
		bool badWindow =
		  ( windowSize < 0 ) &&
		  ( windowSize != -1 );
		if( badWindow ) { parser->setError( "window size" ); }
	}

	// Get the experimental bonus data and options.
	if( !parser->isError() ) {
		experimentalFile = parser->getOptionString( experimentalOptions );
		if( !parser->isError() ) { parser->setOptionDouble( experimentalOffsetOptions, experimentalOffset ); }
		if( !parser->isError() ) { parser->setOptionDouble( experimentalScalingOptions, experimentalScaling ); }
	}
	
	// FD
	// Get the bin size for the histogram.
	if( !parser->isError() ) {
		parser->setOptionDouble( binSizeOptions, binSize );
		if (binSize <= 0.0) {parser->setError( "bin size" );}
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

string Fold::sample_file(string shapefile, int numnuc, int iter) {
	char numstr[21];
	ifstream infile;
	sprintf(numstr, "%d", iter);
	string outname = shapefile + "_boot" + numstr;
	ofstream outfile;
	outfile.open(outname.c_str());
	infile.open(shapefile.c_str());
	double *SHAPEdata = new double [numnuc];
	int index, ridx, i;
	double value;

	for( i=0;i<numnuc;i++ ) { SHAPEdata[i] = -999; }

	while( infile >> index >> value ) { SHAPEdata[index] = value; }

	srand( time(0) );
	for( i=0;i<numnuc;i++ ) {
		ridx = rand() % numnuc;
		if( SHAPEdata[ridx] != -999 ) { SHAPEdata[ridx] += SHAPEdata[ridx]; }
	}

	for( i=0;i<numnuc;i++ ) { outfile << i << " " << SHAPEdata[i] << "\n"; }
	
	outfile.close();
	infile.close();

	return outname;
}
///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void Fold::run() {

	// Create a variable that handles errors.
	int error = 0;
	char numstr[21];
	int b_iter;
	string consFile;

	/*
	 * Use the constructor for RNA that specifies a filename.
	 * Specify type = 2 (sequence file).
	 * isRNA identifies whether the strand is RNA (true) or DNA (false).
	 *
	 * After construction of the strand data structure, create the error checker which monitors for errors.  
	 * Throughout, the error status of the calculation is checked with a variant of the isErrorStatus method, which returns 0 if no error occurred.
	 * The calculation proceeds as long as error = 0.
	 */
	cout << "Initializing nucleic acids..." << flush;
	RNA* strand = new RNA( seqFile.c_str(), 2, isRNA );
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>( strand );
	error = checker->isErrorStatus();
	if( error == 0 ) { cout << "done." << endl; }
	
	/*
	 * FD
 	 * Generate reactivity histogram from training set
 	 * To be finished
     */
	/*if( error == 0 && binSize >= 0 ){
		cout << "Generating reactivity histogram..." << flush;
		
		//int tempError = strand->GenerateHistogram()
		
		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}*/
    

	/*
	 * Set the window size, based on the length of the sequence given as input.
	 * Only do this if window size hasn't been set on the command line.
	 * Use method GetSequenceLength to get the length of the sequence.
	 * The window sizes in relation to the length are hardcoded values.
	 */
	if( windowSize == -1 && error == 0 ) {
		int length = strand->GetSequenceLength();
		windowSize =
			( length > 1200 ) ? 20 :
			( length > 800 ) ? 15 :
			( length > 500 ) ? 11 :
			( length > 300 ) ? 7 :
			( length > 120 ) ? 5 :
			( length > 50 ) ? 3 :
			2;
	}

	/*
	 * Set the temperature using the SetTemperature method.
	 * Only set the temperature if a given temperature doesn't equal the default.
	 * If the temperature does need to be set, use the error checker's isErrorStatus method to check for errors.
	 */
	if( ( error == 0 ) && ( temperature != 310.15 ) ) {

		// Show a message saying that the temperature is being set.
		cout << "Setting temperature..." << flush;

		// Set the temperature and check for errors.
		int tempError = strand->SetTemperature( temperature );
		error = checker->isErrorStatus( tempError );

		// If no error occurred, print a message saying that temperature is set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	/*
	 * Set maximum pairing distance using the ForceMaximumPairingDistance method.
	 */
	if( error == 0 && maxDistance != -1 ) {

		// Show a message saying that the maximum pairing distance is being set.
		cout << "Setting maximum distance between paired nucleotides..." << flush;

		// Set the maximum pairing distance and check for errors.
		int distError = strand->ForceMaximumPairingDistance( maxDistance );
		error = checker->isErrorStatus( distError );

		// If no error occurred, print a message saying that maximum pairing distance was set.
		if( error == 0 ) { cout << "done." << endl; }
	}

	// FD
	/*
	 * set two state option and smooth option
	*/
	if(error == 0) {
		strand->setStateType(twoStateVersion);
		strand->setSmoothVersion(smoothVersion);
		if(error != 0) {
			cout<<"error setting state option and smooth option.\n";
			exit(0);
		}
	}
	
	/*
	 * Add constraints if files have been given for their inclusion.
	 * For folding constraints, use the ReadConstraints method.
	 * For SHAPE constraints, use the ReadSHAPE method.
	 * For single strand offset, use the ReadSSO method.
	 * For double strand offset, use the ReadDSO method.
	 * For experimental pair bonuses, use the ReadExperimentalPairBonus method.
	 * After each constraint type, check the error status of the strand as above.
	 */

	// Determine if constraints should be applied.
	bool applyConstraints =
		( constraintFile != "" ) ||
		( SHAPEFile != "" ) ||
		( DSHAPEFile != "" ) ||
		( CMCTFile != "" ) ||
		( DMSFile != "" ) ||
		( singleOffsetFile != "" ) ||
		( doubleOffsetFile != "" ) ||
		( experimentalFile != "" );

	for( b_iter=0; b_iter <= bootstrap; b_iter++) {
		// If constraints should be applied, do so.
		if( error == 0 && applyConstraints ) {

			// Show a message saying that constraints are being applied.
			cout << "Applying constraints..." << flush;
			int constraintError = 0;

			// Read folding constraints, if applicable.
			if( constraintFile != "" ) {
				constraintError = strand->ReadConstraints( constraintFile.c_str() );
				error = checker->isErrorStatus( constraintError );
			}

			// Read SHAPE constraints
			if( error == 0 && SHAPEFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(SHAPEFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = SHAPEFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, "SHAPE" );
				error = checker->isErrorStatus( constraintError );
			}

			// Read differential SHAPE constraints
			if( error == 0 && DSHAPEFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(DSHAPEFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = DSHAPEFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), Dslope, 0, 0, 0, "diffSHAPE" );
				error = checker->isErrorStatus( constraintError );
			}

			// Read DMS constraints.
			if( error == 0 && DMSFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(DMSFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = DMSFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, "DMS" );
				error = checker->isErrorStatus( constraintError );
			}

			// Read CMCT constraints.
			if( error == 0 && CMCTFile != "" ) {
				if( b_iter > 0 ) { consFile = sample_file(CMCTFile, strand->GetSequenceLength(), b_iter); }
				else { consFile = CMCTFile; }
				constraintError = strand->ReadSHAPE( consFile.c_str(), slope, intercept, slopeSingle, interceptSingle, "CMCT" );
				error = checker->isErrorStatus( constraintError );
			}


			// Read single strand offset, if applicable.
			if( error == 0 && singleOffsetFile != "" ) {
				constraintError = strand->ReadSSO( singleOffsetFile.c_str() );
				error = checker->isErrorStatus( constraintError );
			}

			// Read double strand offset, if applicable.
			if( error == 0 && doubleOffsetFile != "" ) {
				constraintError = strand->ReadDSO( doubleOffsetFile.c_str() );
				error = checker->isErrorStatus( constraintError );
			}

			// Read experimental pair bonus constraints, if applicable.
			if( error == 0 && experimentalFile != "" ) {
				constraintError = strand->ReadExperimentalPairBonus( experimentalFile.c_str(), experimentalOffset, experimentalScaling );
				error = checker->isErrorStatus( constraintError );
			}

			// If no error occurred, print a message saying that constraints were applied.
			if( error == 0 ) { cout << "done." << endl; }
		}

		//Make sure the user isn't using -mfe and -s, these are incompatible.

		if (quickfold&&saveFile!="") {

			error = 1;
			cerr << "Fold stopped.  The -s and -mfe commands are incompatible.\n";

		}

		/*
		 * Fold the single strand using the FoldSingleStrand method.
		 * During calculation, monitor progress using the TProgressDialog class and the Start/StopProgress methods of the RNA class.
		 * Neither of these methods require any error checking.
		 * After the main calculation is complete, use the error checker's isErrorStatus method to check for errors.
		 */
		if( error == 0 ) {

			// Show a message saying that the main calculation has started.
			cout << "Folding single strand..." << flush;

			// Create the progress monitor.
			TProgressDialog* progress = new TProgressDialog();
			strand->SetProgress( *progress );

			// Do the main calculation and check for errors.
			int mainCalcError = strand->FoldSingleStrand( percent, maxStructures, windowSize, saveFile.c_str(), maxLoop, quickfold );
			error = checker->isErrorStatus( mainCalcError );

			// Delete the progress monitor.
			strand->StopProgress();
			delete progress;

			// If no error occurred, print a message saying that the main calculation is done.
			if( error == 0 ) { cout << "done." << endl; }
		}

		/*
		 * Write a CT output file using the WriteCt method.
		 * After writing is complete, use the error checker's isErrorStatus method to check for errors.
		 */
		if( error == 0 ) {

			// Show a message saying that the CT file is being written.
			cout << "Writing output ct file..." << flush;

			// Write the CT file and check for errors.
			sprintf(numstr, "%d", b_iter);
			int writeError;

			if( b_iter > 0 ){ writeError = strand->WriteCt( (ctFile + "_boot" + numstr).c_str() ); }
			
			else { writeError = strand->WriteCt( ctFile.c_str() ); }

			error = checker->isErrorStatus( writeError );

			// If no errors occurred, show a CT file writing completion message.
			if( error == 0 ) { cout << "done." << endl; }
		}
	}

	// Delete the error checker and data structure.
	delete checker;
	delete strand;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	Fold* runner = new Fold();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
