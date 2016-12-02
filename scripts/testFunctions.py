from readMRC import ReadMRC


test = ReadMRC('test-stack.mrc')
#test.read()
#a = test.bfactorReal(0, 150)
print(test.create_coeff_matrix(4))