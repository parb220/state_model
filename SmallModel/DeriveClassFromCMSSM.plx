#/usr/bin/perl

use warnings; 
use strict; 

use Getopt::Long; 
use File::Spec::Functions; 
use File::Path; 

my ($var_file_name, $class_dir); 
my $class_name_suffix = "Derived";
my $help=""; 
# parse command line to get the 
# 	file name containing the variable names
# 	directory to store derived classes
#	suffix of the derived classes
GetOptions("file|f=s" =>\$var_file_name,	# file storing variable names
	   "directory|d=s" => \$class_dir, # directory to store derived classes
	   "name|n=s"=> \$class_name_suffix, # name suffix of the derived classes
	   "help|h" => \$help); 
die "usage: perl DeriveClass.plx -f file -d directory -n name -h\n" if ( $help or not defined $var_file_name or not defined $class_dir );

# read out file containing the variable names
open VAR_FILE, $var_file_name or die "Cannot open $var_file_name for reading the variable names\n"; 
my @file_readout = <VAR_FILE>; 
close VAR_FILE; 

# get the variable names
my @variable; 
foreach my $line (@file_readout)
{
	chomp $line; 
	push(@variable, $line) if $line =~ /^[a-zA-Z]/; 
}

# #debug
# print join("\n", sort(@variable)); 

# Write the derived class 
mkdir ($class_dir) or die $! unless -d $class_dir; 
my ($re_header_file, $re_file, $me_header_file, $me_file, $class_hpp_file, $class_cpp_file, $class_name, $class_name_uppercase); 

##	(1). Rational Expectation .hpp 
$class_name = "RationalExpectationFunction_" . $class_name_suffix; 
$class_hpp_file = $class_name . ".hpp"; 
$re_header_file = catfile($class_dir, $class_hpp_file);
 
open RE_FILE, ">$re_header_file" or die "Cannot open $re_header_file for writing\n";

$class_name_uppercase = $class_name; 
$class_name_uppercase =~ tr/a-z/A-Z/; 
print RE_FILE "#ifndef _", $class_name_uppercase, "_\n"; 
print RE_FILE "#define _", $class_name_uppercase, "_\n"; 

print RE_FILE "\n#include \"CMSSM.hpp\"\n"; 
print RE_FILE "\nclass $class_name : public RationalExpectationFunction\n"; 
print RE_FILE "{\n"; 

print RE_FILE "protected:\n"; 
print RE_FILE "\tvirtual void ConvertXtoParameter(const TDenseVector &x); \n"; 

# variables
foreach (@variable)
{
	print RE_FILE "\tdouble ", $_, ";\n";
}

print RE_FILE "\npublic:\n"; 
print RE_FILE "\tvirtual int convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nE, size_t nExpectation);\n"; 

print RE_FILE "\t$class_name() : RationalExpectationFunction() {}\n"; 
print RE_FILE "\t$class_name(const TDenseVector &p) : RationalExpectationFunction(p) {}\n"; 
print RE_FILE "\t~$class_name() {}\n"; 
print RE_FILE "};\n\n";

print RE_FILE "#endif\n";
close(RE_FILE);

## 	(2) RationalExpectationFunction::ConvertXtoParameter().cpp 
$class_name = "RationalExpectationFunction_" . $class_name_suffix;
$class_cpp_file = $class_name . "_ConvertXtoParameter.cpp"; 
$re_file = catfile($class_dir, $class_cpp_file); 
open RE_FILE, ">$re_file" or die "Cannot open $re_file for writing.\n"; 

$class_hpp_file = $class_name . ".hpp";
print RE_FILE "#include \"$class_hpp_file\"\n"; 
print RE_FILE "\nusing namespace std;\n"; 
$class_name = "RationalExpectationFunction_" . $class_name_suffix;
print RE_FILE "void $class_name :: ConvertXtoParameter(const TDenseVector &x);\n"; 
print RE_FILE "{\n"; 
print RE_FILE "}\n";

close(RE_FILE); 

## 	(3) MeasurementEquationFunction .hpp
$class_name = "MeasurementEquationFunction_" . $class_name_suffix; 
$class_hpp_file = $class_name . ".hpp"; 
$me_header_file = catfile($class_dir, $class_hpp_file); 
open ME_FILE, ">$me_header_file" or die "Cannot open $me_header_file for writing.\n"; 

$class_name_uppercase = $class_name;
$class_name_uppercase =~ tr/a-z/A-Z/;
print ME_FILE "#ifndef _", $class_name_uppercase, "_\n";
print ME_FILE "#define _", $class_name_uppercase, "_\n";

print ME_FILE "\n#include \"CMSSM.hpp\"\n";
print ME_FILE "\nclass $class_name : public MeasurementEquationFunction\n"; 
print ME_FILE "{\n"; 
print ME_FILE "protected:\n"; 
print ME_FILE "\tvirtual void ConvertXtoParameter(const TDenseVector &x)\n"; 
foreach (@variable)
{
	print ME_FILE "\tdouble ", $_, ";\n"; 
}

print ME_FILE "public:\n";
print ME_FILE "\tvirtual int convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nNu);\n"; 
print ME_FILE "\t$class_name() : MeasurementEquationFunction() {}\n"; 
print ME_FILE "\t$class_name(const TDenseVector &p) : MeasurementEquationFunction(p) {}\n"; 
print ME_FILE "\tvirtual ~$class_name() {}\n";
print ME_FILE "};\n\n";

print ME_FILE "#endif\n";
close(ME_FILE);

##	(4) MeasurementEquationFunction::ConvertXtoParameter().cpp
$class_name = "MeasurementEquationFunction_" . $class_name_suffix;
$class_cpp_file = $class_name . "_ConvertXtoParameter.cpp";
$me_file = catfile($class_dir, $class_cpp_file);
open ME_FILE, ">$me_file" or die "Cannot open $me_file for writing.\n";

$class_hpp_file = $class_name . ".hpp";
print ME_FILE "#include \"$class_hpp_file\"\n";
print ME_FILE "\nusing namespace std;\n";
$class_name = "MeasurementEquationFunction_" . $class_name_suffix;
print ME_FILE "void $class_name :: ConvertXtoParameter(const TDenseVector &x)\n";
print ME_FILE "{\n";
print ME_FILE "}\n";

close(RE_FILE);


