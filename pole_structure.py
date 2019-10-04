##################################
##########Import Libraries########
##################################
#Begin
from Label_self import *
import pickle
from copy import deepcopy
#End
###################################

###################################
######Is transformation valid?#####
###################################
def valid(g_freq,p):
  count_p=0
  for i in range(0,len(g_freq)):
    if g_freq[i][p-1]!=0:
      count_p=count_p+1
  #print count_p
  return count_p
#End
###################################

########################################################
##Comarison of two arrays considering the overall sign##
########################################################
#Begin      
def arr_comp_new(r1, r2):
  ovl_sgn = 1
  r2_new = [None]*len(r2)
  val = False
  if two_arr_eq(r1,r2):
    val = True
    ovl_sgn = 1
  for i in range(0,len(r2)):
    r2_new[i] = -r2[i] 
  if two_arr_eq(r1,r2_new):
    val = True
    ovl_sgn = ovl_sgn*(-1)
  
  return val,ovl_sgn 
#End
########################################################  


M=6

#diag_num=14
pol=[]
#diags=[9,12,14,16,19,26,28,29,32,33,34,35]
diags=[111, 141, 142, 181, 198, 199, 204, 208, 210, 246, 266, 268, 306, 349, 353, 379, 383, 413, 444, 464, 489, 496, 502, 509, 538, 539, 552, 567, 588, 612, 628, 629, 630, 633, 637, 650, 668, 669, 670, 671, 676, 683, 693, 694, 701, 714, 724, 763, 770, 777, 786, 787, 797, 817, 828, 831, 837, 858, 861, 862, 865, 869, 875, 884, 885, 894, 895, 896, 912, 943, 946, 961, 973, 986, 993, 997, 1004, 1025, 1033, 1034, 1041, 1057, 1061, 1062, 1067, 1068, 1074, 1075, 1084, 1091, 1092, 1097, 1107, 1111, 1125, 1129, 1130, 1134, 1140, 1146, 1155, 1162, 1163, 1169, 1175, 1176, 1187, 1190, 1194, 1195, 1212, 1219, 1224, 1232, 1239, 1245, 1258, 1268, 1272, 1276, 1278, 1283, 1286, 1288, 1293, 1294, 1296, 1305, 1309, 1311, 1319, 1320, 1321, 1322, 1327, 1329, 1331, 1335, 1337, 1353, 1359, 1365, 1371, 1372, 1373, 1375, 1378, 1381, 1384, 1388, 1390, 1392, 1393, 1401, 1404, 1405, 1409, 1415, 1424, 1425, 1426, 1429, 1430, 1432, 1433, 1435, 1436, 1439, 1440, 1442, 1443, 1445, 1448, 1453, 1454, 1456, 1459, 1461, 1465, 1467, 1468, 1479, 1482, 1484, 1490, 1496, 1497, 1501, 1504, 1506, 1507, 1509, 1514, 1515, 1518, 1519, 1520, 1521, 1522, 1525, 1527, 1529, 1531, 1534, 1539, 1540, 1541, 1543, 1545, 1556, 1559, 1562, 1563, 1564, 1566, 1571, 1572, 1573, 1574, 1576, 1579, 1580, 1583, 1585, 1587, 1592, 1594, 1598, 1599, 1600, 1603, 1606, 1609, 1610, 1616, 1617, 1618, 1620, 1623, 1624, 1629, 1631, 1632, 1636, 1637, 1639, 1640, 1643, 1645, 1646, 1647, 1648, 1649, 1654, 1661, 1662, 1663, 1668, 1675, 1676, 1677, 1686, 1689, 1690, 1694, 1701, 1706, 1707, 1711, 1714, 1717, 1719, 1722, 1723, 1724, 1725, 1726, 1727, 1728, 1729, 1732, 1734, 1735, 1736, 1737, 1740, 1742, 1745, 1747, 1750, 1754, 1755, 1757, 1759, 1760, 1764, 1765, 1767, 1769, 1770, 1771, 1774, 1775, 1776, 1777, 1779, 1782, 1783, 1789, 1790, 1791, 1793, 1794, 1796, 1798, 1800, 1801, 1802, 1803, 1805, 1806, 1807, 1809, 1810, 1812, 1813, 1814, 1815, 1817, 1818, 1819, 1820, 1822, 1823, 1824, 1826, 1828, 1829, 1832, 1833, 1836, 1837, 1838, 1839, 1840, 1842, 1843, 1844, 1845, 1846, 1847, 1848, 1849, 1850, 1851, 1852, 1853, 1854, 1855, 1856, 1858, 1861, 1862, 1863, 1864, 1865, 1866, 1867, 1868, 1869, 1870, 1871, 1872, 1873, 1874, 1875, 1877, 1881, 1882, 1883, 1885, 1886, 1887, 1888, 1889, 1890, 1891, 1893, 1894, 1896]
for i in diags:
  file_name='m_'+str(M)+'_num_'+str(i)+'.graphml'
  G=load(file_name)
  #print G
  #print G.es["F_or_B"]
  #print G.es["INT_or_EXT"]

  pol_str_i=[]
  pol_str=[]
  #ext_max=0
  for j in range (0,100000):
    reset_g(G,M)
    label_abs_ran(G,M,20,20)
    labels=deepcopy(G.es["Label"])
    ami_in = AMI_Input(labels,M)
    #print 'ami_in = ', ami_in
    pol_str_i.append(valid(ami_in,0))
    pol_str_i.append(valid(ami_in,1))
    pol_str_i.append(valid(ami_in,2))
    pol_str_i.append(valid(ami_in,3))
    pol_str_i.append(valid(ami_in,4))
    pol_str_i.append(valid(ami_in,5))
    pol_str_i.sort()
    pol_str_i.append(valid(ami_in,6))
    if pol_str_i not in pol_str:
      pol_str.append(pol_str_i)
  
    pol_str_i=[]
  pol.append(pol_str)
  print 'graph --> ', i
  file_name_str = 'pole_str_'+str(i)
  with open('results_6/'+file_name_str,"wb") as fp:
    pickle.dump(pol_str,fp)
  #print pol_str
  file_name_len_str = 'len_pole_str_'+str(i)
  with open('results_6/'+file_name_len_str,"wb") as fp:
    pickle.dump(len(pol_str),fp)
  print len(pol_str)
  #with open('results_6/'+file_name_str,"rb") as fp:
    #label=pickle.load(fp)
  #print 'labels = ', label
  #with open('results_6/'+file_name_len_str,"rb") as fp:
    #label=pickle.load(fp)
  #print 'labels = ', label
#print 'pol = ', pol
'''
for j in range (0,2):
  for k in range (j+1,2):
    count=0
    for i in range(0,len(pol[j])):
      if pol[j][i] in pol[k]:
        count=count+1
    print j+1, k+1
    print 'count = ', count
    
#count=0 
#for i in range(0,len(pol[0])):
  #if pol[0][i] in pol[3]:
    #count=count+1
#print 'count = ', count 
'''
