# Input data with the number of mortalities per day, the herd size and the time of cull
# Output as lists

remtimes1 <- c(1,4,7,8,12,13,rep(15,4),rep(16,6),rep(17,5),rep(18,6),19,rep(20,8),22,rep(23,11),rep(24,13))
M1 <- 1614
Tfinal1 <- 24.01
remtimes2 <- c(2,3,4,5,5,6,6,6,7,rep(8,9),rep(9,6),rep(10,4),rep(11,4),rep(12,11),rep(13,3),rep(14,5),rep(15,7),rep(16,2),rep(17,3))
Tfinal2 <- 17.01
M2 <- 1949
remtimes3 <- c(5,6,7,8,8,9,9,rep(10,5),rep(11,4),rep(12,9),13,rep(14,3),rep(15,3),rep(16,5),rep(17,6))
Tfinal3 <- 17.01
M3 <- 1753
remtimes4 <- c(1,2,2,rep(3,5),rep(4,7),rep(5,5),rep(6,7),rep(7,7),rep(8,8),rep(9,4),rep(10,5),rep(11,4),rep(12,10),rep(13,2),rep(14,14),rep(15,18),rep(16,3))
Tfinal4 <- 16.01
M4 <- 1833
remtimes5 <- c(rep(3,4),rep(4,6),rep(5,4),rep(6,4),rep(7,6),rep(8,8),rep(9,6),10,rep(11,3),rep(12,5),rep(13,9),rep(14,8),rep(15,19))
Tfinal5 <- 15.01
M5 <- 1320
remtimes6 <- c(5,6,7,7,8,8,9,12,13,rep(14,8))
Tfinal6 <- 14.01
M6 <- 600
remtimes7 <- c(1,1,3,5,5,6,7,7,7,7,8,8,rep(9,8),rep(10,12),rep(11,8),rep(12,6),rep(13,5))
Tfinal7 <- 14.01
M7 <- 600
remtimes8 <- c(6,6,8,8,8,9,9,9,10,10,10,10,11,11,11,13)
Tfinal8 <- 14.01
M8 <- 600
remtimes9 <- c(rep(2,7),rep(3,7),rep(4,3),rep(5,2),rep(7,5),8,rep(10,6),12,rep(13,4),rep(15,10),rep(16,36),rep(17,24),rep(18,17),rep(19,18),rep(20,7),rep(21,40),rep(22,42),rep(23,31),rep(24,12))
Tfinal9 <- 24.01
M9 <- 2145

herd1 <- list(remtimes=remtimes1, M=M1, Tfinal=Tfinal1)
herd2 <- list(remtimes=remtimes2, M=M2, Tfinal=Tfinal2)
herd3 <- list(remtimes=remtimes3, M=M3, Tfinal=Tfinal3)
herd4 <- list(remtimes=remtimes4, M=M4, Tfinal=Tfinal4)
herd5 <- list(remtimes=remtimes5, M=M5, Tfinal=Tfinal5)
herd6 <- list(remtimes=remtimes6, M=M6, Tfinal=Tfinal6)
herd7 <- list(remtimes=remtimes7, M=M7, Tfinal=Tfinal7)
herd8 <- list(remtimes=remtimes8, M=M8, Tfinal=Tfinal8)
herd9 <- list(remtimes=remtimes9, M=M9, Tfinal=Tfinal9)