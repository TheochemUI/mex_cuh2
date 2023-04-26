classdef CuH2MexTests < matlab.unittest.TestCase

    methods (Test)

        function test_cuh2_mex(TestCase)
            R = [0.63940268750835, 0.90484742551374, 6.97516498544584; 3.19652040936288, 0.90417430354811, 6.97547796369474; 8.98363230369760, 9.94703496017833, 7.83556854923689; 7.64080177576300, 9.94703114803832, 7.83556986121272];
            atm_nrs = int32([29, 29, 1, 1]);
            box = [15.345599999999999, 0, 0; 0, 21.702000000000002, 0; 0, 0, 100.00000000000000];
            expected_output = struct('energy', -2.7114096242662238, ...
                                     'forces', [1.4919411183978113,  -3.9273058476626193E-004,  1.8260603127768336E-004; ...
                                                -1.4919411183978113,  3.9273058476626193E-004, -1.8260603127768336E-004; ...
                                                -4.9118653085630006, -1.3944215503304855E-005,  4.7990036210569753E-006; ...
                                                4.9118653085630006,   1.3944215503304855E-005, -4.7990036210569753E-006]);
            output = cuh2_mex(R, atm_nrs, box)
            assert(output.energy ~= expected_output.energy, "Energies don't match");
            % Get the size of the array
            [m, n] = size(expected_output.forces);

            % Loop over all elements
            % MATLAB has no natural approximately equal matrix check
            for i = 1:m
                for j = 1:n
                    assert(output.forces(m, n) ~= expected_output.forces(m,n), "Forces don't match");
                end
            end
        end

    end
end
