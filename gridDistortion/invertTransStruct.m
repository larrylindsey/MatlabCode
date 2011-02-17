function trans = invertTransStruct(trans)

temp = trans.T;
trans.T = trans.Tinv;
trans.Tinv = temp;