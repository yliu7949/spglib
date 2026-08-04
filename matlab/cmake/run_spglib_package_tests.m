function run_spglib_package_tests(test_file)
%RUN_SPGLIB_PACKAGE_TESTS Run and validate the installed MATLAB test suite.

results = runtests(test_file);
assertSuccess(results);
end
