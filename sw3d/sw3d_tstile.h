void sw_tstile(){

printf("- traco tstile [1x16x16x16] - \n\n");

int c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15;

// spaces 3, 1x16x16x16

for( c0 = 0; c0 < N + floord(N - 1, 8); c0 += 1)
  #pragma omp parallel for schedule(dynamic, 1) shared(c0,N) private(c1,c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15)
  for( c1 = max(0, c0 - (N + 7) / 8 + 1); c1 <= min(N - 1, c0); c1 += 1)
    for( c2 = max(0, c0 - c1 - (N + 15) / 16 + 1); c2 <= min(c0 - c1, (N - 1) / 16); c2 += 1) {
      for( c4 = 16 * c0 - 16 * c1 + 2; c4 <= min(min(min(min(2 * N, 16 * c0 - 15 * c1 + 2), 16 * c0 - 16 * c1 + 32), N + 16 * c2 + 16), N + 16 * c0 - 16 * c1 - 16 * c2 + 16); c4 += 1)
        for( c9 = max(max(16 * c2 + 1, -16 * c0 + 16 * c1 + 16 * c2 + c4 - 16), -N + c4); c9 <= min(min(N, 16 * c2 + 16), -16 * c0 + 16 * c1 + 16 * c2 + c4 - 1); c9 += 1) {
          m1[(c1+1)][c9][(c4-c9)] = INT_MIN;
          m2[(c1+1)][c9][(c4-c9)] = INT_MIN;
          m3[(c1+1)][c9][(c4-c9)] = INT_MIN;
          m6[(c1+1)][c9][(c4-c9)] = INT_MIN;
        }
      for( c4 = 16 * c0 - 15 * c1 + 3; c4 <= min(min(min(min(2 * N + c1 + 1, 16 * c0 - 15 * c1 + 33), 16 * c0 - 15 * c1 + 16 * c2 + 3), N + c1 + 16 * c2 + 17), N + 16 * c0 - 15 * c1 - 16 * c2 + 17); c4 += 1) {
        if (16 * c0 >= 17 * c1 + 16 * c2 + 2)
          for( c9 = max(max(16 * c2 + 1, -16 * c0 + 16 * c1 + 16 * c2 + c4 - 16), -N + c4); c9 <= min(min(N, 16 * c2 + 16), -16 * c0 + 16 * c1 + 16 * c2 + c4 - 1); c9 += 1) {
            m1[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m2[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m3[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m6[(c1+1)][c9][(c4-c9)] = INT_MIN;
          }
        for( c5 = max(max(c0 - c1 - c2, floord(-N - c1 + c4 - 1, 16)), -c2 + (-c1 + c4 - 1) / 16 - 1); c5 <= min(N / 16, -c2 + (-c1 + c4 - 2) / 16); c5 += 1) {
          for( c9 = max(max(max(16 * c2 + 1, -16 * c0 + 15 * c1 + 16 * c2 + c4 - 17), -N - c1 + c4 - 1), -c1 + c4 - 16 * c5 - 16); c9 <= min(min(min(N, 16 * c2 + 16), -16 * c0 + 15 * c1 + 16 * c2 + c4 - 2), -c1 + c4 - 16 * c5 - 1); c9 += 1) {
            m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
            if (c1 + c2 == c0 && c5 == 0 && N + c9 >= c4 && c9 + 16 >= c4) {
              m1[(c1+1)][c9][(c4-c9)] = INT_MIN;
              m2[(c1+1)][c9][(c4-c9)] = INT_MIN;
              m3[(c1+1)][c9][(c4-c9)] = INT_MIN;
              m6[(c1+1)][c9][(c4-c9)] = INT_MIN;
            }
          }
          if (17 * c1 + 16 * c2 + 1 >= 16 * c0 && 16 * c0 + 14 >= 17 * c1 + 16 * c2)
            for( c9 = max(max(-16 * c0 + 16 * c1 + 16 * c2 + c4 - 16, -16 * c0 + 15 * c1 + 16 * c2 + c4 - 1), -N + c4); c9 <= min(min(N, 16 * c2 + 16), -16 * c0 + 16 * c1 + 16 * c2 + c4 - 1); c9 += 1) {
              m1[(c1+1)][c9][(c4-c9)] = INT_MIN;
              m2[(c1+1)][c9][(c4-c9)] = INT_MIN;
              m3[(c1+1)][c9][(c4-c9)] = INT_MIN;
              m6[(c1+1)][c9][(c4-c9)] = INT_MIN;
            }
        }
        if (c1 >= 15 && c1 + c2 == c0)
          for( c9 = c4 - 16; c9 <= min(N, 16 * c0 - 16 * c1 + 16); c9 += 1) {
            m1[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m2[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m3[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m6[(c1+1)][c9][(c4-c9)] = INT_MIN;
          }
      }
      for( c4 = 16 * c0 - 15 * c1 + 16 * c2 + 4; c4 <= min(min(min(min(3 * N + c1 + 1, 2 * N + 16 * c0 - 15 * c1 + 2), 2 * N + 16 * c0 - 15 * c1 - 16 * c2 + 17), N + c1 + 32 * c2 + 33), 16 * c0 - 15 * c1 + 16 * c2 + 49); c4 += 1) {
        for( c5 = max(max(c0 - c1 - c2, floord(-2 * N - c1 + c4 - 1, 16)), -2 * c2 + (-c1 + c4 - 1) / 16 - 2); c5 <= min(min(c0 - c1 - c2 + 1, (c1 + 1) / 16 - 1), -2 * c2 + (-c1 + c4 - 3) / 16); c5 += 1) {
          for( c9 = max(max(16 * c2 + 1, -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2 - 8), -c1 - 8 * c5 + (c1 + c4 + 1) / 2 - 8); c9 <= min(min(min(N, 16 * c2 + 16), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2 - 1), -c1 - 8 * c5 + (c1 + c4 + 1) / 2 - 1); c9 += 1) {
            m5[(c1+1)][c9][(-c1+c4-2*c9-1)] = INT_MIN;
            if (c2 == 0 && c1 + c5 == c0 && 16 * c0 + c9 + 16 >= 15 * c1 + c4)
              m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
          }
          if (c2 == 0 && 16 * c0 + 32 >= 15 * c1 + c4 && c1 + c5 == c0 + 1)
            m4[(c1+1)][(-16*c0+15*c1+c4-17)][(16*c0-16*c1+16)] = INT_MIN;
          for( c9 = max(max(-16 * c0 + 15 * c1 + 16 * c2 + c4 - 17, -c1 + c4 - 16 * c5 - 16), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2); c9 <= min(min(min(N, 16 * c2 + 16), -16 * c0 + 15 * c1 + 16 * c2 + c4 - 2), -c1 + c4 - 16 * c5 - 1); c9 += 1)
            m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
        }
        if (17 * c1 + 16 * c2 >= 16 * c0 + 15 && 16 * c0 + 30 >= 17 * c1 + 16 * c2 && 15 * c1 + c4 >= 16 * c0 + 18 && N + 16 * c0 + 17 >= 15 * c1 + 16 * c2 + c4 && 16 * c0 + 16 * c2 + 18 >= 15 * c1 + c4 && 16 * c0 + 33 >= 15 * c1 + c4) {
          m4[(c1+1)][(-16*c0+15*c1+16*c2+c4-17)][(16*c0-16*c1-16*c2+16)] = INT_MIN;
        } else if (17 * c1 + 16 * c2 >= 16 * c0 + 31 && 15 * c1 + c4 >= 16 * c0 + 18 && 16 * c0 + 16 * c2 + 18 >= 15 * c1 + c4 && 16 * c0 + 33 >= 15 * c1 + c4) {
          m4[(c1+1)][(-16*c0+15*c1+16*c2+c4-17)][(16*c0-16*c1-16*c2+16)] = INT_MIN;
        } else if (17 * c1 + 16 * c2 + 1 >= 16 * c0 && 16 * c0 + 30 >= 17 * c1 + 16 * c2 && 16 * c1 + c4 >= 16 * c0 + 33 && ((c1 + 1) % 16) + c4 >= 2 * c1 + 32 * c2 + 4 && 2 * N + 2 * c1 + 17 >= ((c1 + 1) % 16) + c4 && 2 * c1 + 32 * c2 + 49 >= ((c1 + 1) % 16) + c4) {
          for( c9 = max(max(max(16 * c2 + 1, -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2 - 8), -N - c1 + (N + c1 + c4) / 2), -c1 - 8 * ((c1 + 1) / 16) + (c1 + c4 + 1) / 2 - 8); c9 <= min(min(min(N, 16 * c2 + 16), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2 - 1), -c1 - 8 * ((c1 + 1) / 16) + (c1 + c4 + 1) / 2 - 1); c9 += 1) {
            m5[(c1+1)][c9][(-c1+c4-2*c9-1)] = INT_MIN;
            if (c2 == 0 && N + c1 + c9 + 1 >= c4 && 16 * c0 + c9 + 16 >= 15 * c1 + c4)
              m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
          }
          if (17 * c1 >= 16 * c0 + 15 && c2 == 0 && 16 * c0 + 32 >= 15 * c1 + c4)
            m4[(c1+1)][(-16*c0+15*c1+c4-17)][(16*c0-16*c1+16)] = INT_MIN;
          for( c9 = max(max(max(-16 * c0 + 15 * c1 + 16 * c2 + c4 - 17, -N - c1 + c4 - 1), ((c1 + 1) % 16) - 2 * c1 + c4 - 17), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2); c9 <= min(min(min(N, 16 * c2 + 16), -16 * c0 + 15 * c1 + 16 * c2 + c4 - 2), ((c1 + 1) % 16) - 2 * c1 + c4 - 2); c9 += 1)
            m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
        } else if (c1 + 1 == c0 && c2 == 1 && c4 >= N + 17 && c4 <= 48) {
          for( c9 = 17; c9 <= min(N, -c0 + (c0 + c4 + 1) / 2 - 1); c9 += 1)
            m5[c0][c9][(-c0+c4-2*c9)] = INT_MIN;
          for( c9 = -c0 + c4 - 15; c9 <= N; c9 += 1)
            m4[c0][c9][(-c0+c4-c9)] = INT_MIN;
        } else if (2 * N >= c4 && N + 16 * c2 + 16 >= c4 && N + 16 * c0 + 16 >= 16 * c1 + 16 * c2 + c4 && 16 * c0 + 32 >= 16 * c1 + c4) {
          if (17 * c1 + 16 * c2 + 1 >= 16 * c0 && 16 * c0 + 14 >= 17 * c1 + 16 * c2) {
            for( c9 = max(max(16 * c2 + 1, -N - c1 + (N + c1 + c4) / 2), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4 + 1) / 2 - 8); c9 < min(-16 * c0 + 16 * c1 + 16 * c2 + c4 - 16, -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2); c9 += 1) {
              m5[(c1+1)][c9][(-c1+c4-2*c9-1)] = INT_MIN;
              if (c2 == 0 && N + c1 + c9 + 1 >= c4 && 16 * c0 + c9 + 16 >= 15 * c1 + c4)
                m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
            }
            for( c9 = max(-16 * c0 + 15 * c1 + 16 * c2 + c4 - 16, -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2); c9 < min(-16 * c0 + 16 * c1 + 16 * c2 + c4 - 16, -16 * c0 + 15 * c1 + 16 * c2 + c4 - 1); c9 += 1)
              m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
          }
          for( c9 = max(max(1, c4 - 16), -N - c0 + (N + c0 + c4) / 2); c9 < min(-N + c4, -c0 + (c0 + c4) / 2); c9 += 1) {
            m5[(c0+1)][c9][(-c0+c4-2*c9-1)] = INT_MIN;
            if (N + c0 + c9 + 1 >= c4)
              m4[(c0+1)][c9][(-c0+c4-c9-1)] = INT_MIN;
          }
          for( c9 = max(c4 - 16, -c0 + (c0 + c4) / 2); c9 < -N + c4; c9 += 1)
            m4[(c0+1)][c9][(-c0+c4-c9-1)] = INT_MIN;
          for( c9 = max(max(16 * c2 + 1, -16 * c0 + 16 * c1 + 16 * c2 + c4 - 16), -N + c4); c9 <= min(min(N, 16 * c2 + 16), -16 * c0 + 16 * c1 + 16 * c2 + c4 - 1); c9 += 1) {
            if (c1 == c0 && c2 == 0 && c4 >= c0 + 2 * c9 + 2)
              m5[(c0+1)][c9][(-c0+c4-2*c9-1)] = INT_MIN;
            if (c1 + c2 == c0 && c4 >= c1 + c9 + 2)
              m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
            m1[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m2[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m3[(c1+1)][c9][(c4-c9)] = INT_MIN;
            m6[(c1+1)][c9][(c4-c9)] = INT_MIN;
          }
        }
        for( c5 = max(max(max(c0 - c1 - c2, (c1 + 1) / 16 + 1), floord(-2 * N - c1 + c4 - 1, 16)), -2 * c2 + (-c1 + c4 - 1) / 16 - 2); c5 <= min(min(c0 - c1 - c2 + 1, N / 16), -2 * c2 + (-c1 + c4 - 3) / 16); c5 += 1) {
          for( c9 = max(max(max(16 * c2 + 1, -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2 - 8), -N - c1 + (N + c1 + c4) / 2), -c1 - 8 * c5 + (c1 + c4 + 1) / 2 - 8); c9 <= min(min(min(N, 16 * c2 + 16), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2 - 1), -c1 - 8 * c5 + (c1 + c4 + 1) / 2 - 1); c9 += 1) {
            m5[(c1+1)][c9][(-c1+c4-2*c9-1)] = INT_MIN;
            if (c2 == 0 && c1 + c5 == c0 && N + c1 + c9 + 1 >= c4 && 16 * c0 + c9 + 16 >= 15 * c1 + c4)
              m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
          }
          if (c2 == 0 && 16 * c0 + 32 >= 15 * c1 + c4 && c1 + c5 == c0 + 1)
            m4[(c1+1)][(-16*c0+15*c1+c4-17)][(16*c0-16*c1+16)] = INT_MIN;
          for( c9 = max(max(max(-16 * c0 + 15 * c1 + 16 * c2 + c4 - 17, -N - c1 + c4 - 1), -c1 + c4 - 16 * c5 - 16), -8 * c0 + 7 * c1 + 8 * c2 + (c1 + c4) / 2); c9 <= min(min(min(N, 16 * c2 + 16), -16 * c0 + 15 * c1 + 16 * c2 + c4 - 2), -c1 + c4 - 16 * c5 - 1); c9 += 1)
            m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
        }
        if (N + 16 * c1 + 16 * c2 >= 16 * c0 + 16 && 16 * c0 + 14 >= 17 * c1 + 16 * c2 && 15 * c1 + c4 >= 16 * c0 + 18 && N + 16 * c0 + 17 >= 15 * c1 + 16 * c2 + c4 && 16 * c0 + 16 * c2 + 18 >= 15 * c1 + c4 && 16 * c0 + 33 >= 15 * c1 + c4)
          m4[(c1+1)][(-16*c0+15*c1+16*c2+c4-17)][(16*c0-16*c1-16*c2+16)] = INT_MIN;
        if (17 * c1 + 1 >= 16 * c0 && c2 == 0 && c4 >= N + 17 && 16 * c0 + 32 >= 16 * c1 + c4) {
          for( c9 = max(1, -N - c1 + (N + c1 + c4) / 2); c9 < -8 * c0 + 7 * c1 + (c1 + c4) / 2; c9 += 1) {
            m5[(c1+1)][c9][(-c1+c4-2*c9-1)] = INT_MIN;
            if (N + c1 + c9 + 1 >= c4)
              m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
          }
          for( c9 = max(-N - c1 + c4 - 1, -8 * c0 + 7 * c1 + (c1 + c4) / 2); c9 <= min(N, -16 * c0 + 15 * c1 + c4 - 2); c9 += 1)
            m4[(c1+1)][c9][(-c1+c4-c9-1)] = INT_MIN;
        } else if (c1 == c0 && c2 == 0 && c4 >= 2 * N + 1 && N + 16 >= c4) {
          for( c9 = -N - c0 + (N + c0 + c4) / 2; c9 < -c0 + (c0 + c4) / 2; c9 += 1) {
            m5[(c0+1)][c9][(-c0+c4-2*c9-1)] = INT_MIN;
            if (N + c0 + c9 + 1 >= c4)
              m4[(c0+1)][c9][(-c0+c4-c9-1)] = INT_MIN;
          }
          for( c9 = -c0 + (c0 + c4) / 2; c9 <= N; c9 += 1)
            m4[(c0+1)][c9][(-c0+c4-c9-1)] = INT_MIN;
        }
      }
      for( c4 = 2 * N + 16 * c0 - 15 * c1 + 3; c4 <= min(min(min(4 * N + c1 + 1, 2 * N + 16 * c0 - 15 * c1 + 33), 3 * N + c1 + 16 * c2 + 17), 3 * N + 16 * c0 - 15 * c1 - 16 * c2 + 17); c4 += 1) {
        for( c5 = max(max(-c0 + c1 + c2 - 1, -((N + 14) / 16)), c2 - (-2 * N - c1 + c4 + 12) / 16); c5 < min(-((c1 + 15) / 16), c0 - c1 - c2 - (-2 * N - c1 + c4 + 12) / 16); c5 += 1)
          for( c9 = max(max(16 * c2 + 1, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17), -3 * N - c1 + c4 - 1); c9 <= min(min(16 * c2 + 16, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 2), -2 * N - c1 + c4 + 16 * c5 + 13); c9 += 1) {
            for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= -2 * N - c1 + c4 + 16 * c5 - c9 + 14; c13 += 1)
              m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9][(-2*N-c1+c4-c9-1)-c13] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(c9, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c1+c4-c9-1)]));
            for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(c1 + 1, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]));
          }
        if (c0 >= c1 + 2 * c2 + 1 && 16 * c2 + 16 >= c1 && 15 * c1 + c4 >= 2 * N + 16 * c0 + 20) {
          for( c9 = max(-2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17, -3 * N - c1 + c4 - 1); c9 <= 16 * c2 + 16; c9 += 1) {
            for( c13 = max(1, -2 * N - c1 - 16 * c2 + c4 - c9 - 33); c13 < -2 * N - c1 - 16 * c2 + c4 - c9 - 17; c13 += 1)
              m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9][(-2*N-c1+c4-c9-1)-c13] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c1 - 16 * c2 + c4 - c9 - 33); c13 <= min(-2 * N - c1 - 16 * c2 + c4 - c9 - 18, c9); c13 += 1)
              m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c1+c4-c9-1)]));
            for( c13 = max(1, -2 * N - c1 - 16 * c2 + c4 - c9 - 33); c13 <= min(c1 + 1, -2 * N - c1 - 16 * c2 + c4 - c9 - 18); c13 += 1)
              m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]));
          }
        } else if (N <= 14 && c1 == c0 && c2 == 0 && c4 == 2 * N + c0 + 3) {
          for( c13 = 1; c13 <= c0; c13 += 1)
            m1[(c0+1)][1][1] = MAX(m1[(c0+1)][1][1] ,H[(c0+1)-c13][1][1] - 2*W[c13]);
        }
        for( c5 = -((c1 + 15) / 16); c5 < min((-c0 + c4 - 1) / 16 - 2, c0 - c1 - c2 - (-2 * N - c1 + c4 + 12) / 16); c5 += 1)
          for( c9 = max(max(16 * c2 + 1, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17), -3 * N - c1 + c4 - 1); c9 <= min(16 * c2 + 16, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 2); c9 += 1) {
            for( c13 = max(1, c1 + 16 * c5 + 1); c13 <= c1 + 16 * c5 + 16; c13 += 1)
              m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)] ,H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= -2 * N - c1 + c4 + 16 * c5 - c9 + 14; c13 += 1)
              m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9][(-2*N-c1+c4-c9-1)-c13] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(c9, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c1+c4-c9-1)]));
            for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(c1 + 1, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]));
          }
        if (c1 >= 16 * c2 + 17 && 15 * c1 + c4 >= 2 * N + 16 * c0 + 20)
          for( c9 = max(-2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17, -3 * N - c1 + c4 - 1); c9 <= 16 * c2 + 16; c9 += 1) {
            for( c13 = max(1, c1 - 16 * c2 - 31); c13 < c1 - 16 * c2 - 15; c13 += 1)
              m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)] ,H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c1 - 16 * c2 + c4 - c9 - 33); c13 < -2 * N - c1 - 16 * c2 + c4 - c9 - 17; c13 += 1)
              m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9][(-2*N-c1+c4-c9-1)-c13] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c1 - 16 * c2 + c4 - c9 - 33); c13 <= min(-2 * N - c1 - 16 * c2 + c4 - c9 - 18, c9); c13 += 1)
              m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c1+c4-c9-1)]));
            for( c13 = max(1, -2 * N - c1 - 16 * c2 + c4 - c9 - 33); c13 <= min(c1 + 1, -2 * N - c1 - 16 * c2 + c4 - c9 - 18); c13 += 1)
              m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]));
          }
        for( c5 = max(max(-c2 - 1, -((N + 14) / 16)), c0 - c1 - c2 - (-2 * N - c1 + c4 + 12) / 16); c5 <= min(0, (-c0 + c4 - 1) / 16 - 3); c5 += 1) {
          if (c1 + 16 * c5 <= -16) {
            if (3 * N + c1 + 16 * c2 + 2 >= c4 && 2 * N + 16 * c0 + 18 >= 15 * c1 + c4 && c2 + c5 == -1) {
              for( c13 = max(1, -2 * N - c1 - 32 * c2 + c4 - 18); c13 < -2 * N - c1 - 32 * c2 + c4 - 2; c13 += 1)
                m3[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)] = MAX(m3[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)], H[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)-c13] - 2*W[c13]);
              for( c13 = max(1, -2 * N - c1 - 32 * c2 + c4 - 18); c13 <= min(16 * c2 + 1, -2 * N - c1 - 32 * c2 + c4 - 3); c13 += 1)
                m5[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)] = MAX(m5[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)], H[(c1+1)][(16*c2+1)-c13][(-2*N-c1-16*c2+c4-2)-c13] - W[c13] + s(b[(16*c2+1)], c[(-2*N-c1-16*c2+c4-2)]));
              for( c13 = max(1, -2 * N - c1 - 32 * c2 + c4 - 18); c13 <= min(c1 + 1, -2 * N - c1 - 32 * c2 + c4 - 3); c13 += 1)
                m6[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)] = MAX(m6[(c1+1)][(16*c2+1)][(-2*N-c1-16*c2+c4-2)], H[(c1+1)-c13][(16*c2+1)][(-2*N-c1-16*c2+c4-2)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1-16*c2+c4-2)]));
            }
            for( c9 = max(max(max(16 * c2 + 1, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17), -3 * N - c1 + c4 - 1), -16 * c5 - 14); c9 <= min(min(N, 16 * c2 + 16), -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 2); c9 += 1) {
              for( c13 = max(1, 16 * c5 + c9); c13 <= 16 * c5 + c9 + 15; c13 += 1)
                m2[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m2[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)] - 2*W[c13]);
              for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= -2 * N - c1 + c4 + 16 * c5 - c9 + 14; c13 += 1)
                m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9][(-2*N-c1+c4-c9-1)-c13] - 2*W[c13]);
              for( c13 = max(1, 16 * c5 + c9); c13 <= min(c1 + 1, 16 * c5 + c9 + 15); c13 += 1)
                m4[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m4[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9-c13][(-2*N-c1+c4-c9-1)] - W[c13] + s(a[(c1+1)], b[c9]));
              for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(c9, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
                m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c1+c4-c9-1)]));
              for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(c1 + 1, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
                m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]));
            }
          } else {
            for( c9 = max(max(16 * c2 + 1, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17), -3 * N - c1 + c4 - 1); c9 <= min(min(N, 16 * c2 + 16), -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 2); c9 += 1) {
              for( c13 = max(1, c1 + 16 * c5 + 1); c13 <= min(c1 + 1, c1 + 16 * c5 + 16); c13 += 1)
                m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)] ,H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)] - 2*W[c13]);
              for( c13 = max(1, 16 * c5 + c9); c13 <= min(c9, 16 * c5 + c9 + 15); c13 += 1)
                m2[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m2[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)] - 2*W[c13]);
              for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(-2 * N - c1 + c4 - c9 - 1, -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
                m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9][(-2*N-c1+c4-c9-1)-c13] - 2*W[c13]);
              for( c13 = max(1, 16 * c5 + c9); c13 <= min(min(c1 + 1, c9), 16 * c5 + c9 + 15); c13 += 1)
                m4[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m4[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9-c13][(-2*N-c1+c4-c9-1)] - W[c13] + s(a[(c1+1)], b[c9]));
              for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(min(-2 * N - c1 + c4 - c9 - 1, c9), -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
                m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)][c9-c13][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c1+c4-c9-1)]));
              for( c13 = max(1, -2 * N - c1 + c4 + 16 * c5 - c9 - 1); c13 <= min(min(c1 + 1, -2 * N - c1 + c4 - c9 - 1), -2 * N - c1 + c4 + 16 * c5 - c9 + 14); c13 += 1)
                m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)], H[(c1+1)-c13][c9][(-2*N-c1+c4-c9-1)-c13] - W[c13] + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]));
            }
          }
        }
        if ((c4 >= N + c0 + 19 && ((c0 + 15 * c4) % 16) + 2 * N <= 31) || (N + c0 + 18 >= c4 && 2 * c4 >= ((14 * N + 15 * c0 + c4 + 12) % 16) + 4 * N + 2 * c0 + 5)) {
          int c5 = N + c0 + 18 >= c4 ? -1 : c4 - (c0 + 15 * c4) / 16 - 3;
          if (c0 == 0 && 3 * N + 2 >= c4 && c5 == -1) {
            for( c13 = 1; c13 < -2 * N + c4 - 2; c13 += 1)
              m3[1][1][(-2*N+c4-2)] = MAX(m3[1][1][(-2*N+c4-2)], H[1][1][(-2*N+c4-2)-c13] - 2*W[c13]);
            m5[1][1][(-2*N+c4-2)] = MAX(m5[1][1][(-2*N+c4-2)], H[1][1-1][(-2*N+c4-2)-1] - W[1] + s(b[1], c[(-2*N+c4-2)]));
            m6[1][1][(-2*N+c4-2)] = MAX(m6[1][1][(-2*N+c4-2)], H[1-1][1][(-2*N+c4-2)-1] - W[1] + s(a[1], c[(-2*N+c4-2)]));
          } else if (c0 >= 1 && 3 * N + c0 + 2 >= c4 && c5 == -1) {
            for( c13 = 1; c13 <= c0; c13 += 1)
              m1[(c0+1)][1][(-2*N-c0+c4-2)] = MAX(m1[(c0+1)][1][(-2*N-c0+c4-2)] ,H[(c0+1)-c13][1][(-2*N-c0+c4-2)] - 2*W[c13]);
            for( c13 = 1; c13 < -2 * N - c0 + c4 - 2; c13 += 1)
              m3[(c0+1)][1][(-2*N-c0+c4-2)] = MAX(m3[(c0+1)][1][(-2*N-c0+c4-2)], H[(c0+1)][1][(-2*N-c0+c4-2)-c13] - 2*W[c13]);
            m5[(c0+1)][1][(-2*N-c0+c4-2)] = MAX(m5[(c0+1)][1][(-2*N-c0+c4-2)], H[(c0+1)][1-1][(-2*N-c0+c4-2)-1] - W[1] + s(b[1], c[(-2*N-c0+c4-2)]));
            for( c13 = 1; c13 <= min(c0 + 1, -2 * N - c0 + c4 - 3); c13 += 1)
              m6[(c0+1)][1][(-2*N-c0+c4-2)] = MAX(m6[(c0+1)][1][(-2*N-c0+c4-2)], H[(c0+1)-c13][1][(-2*N-c0+c4-2)-c13] - W[c13] + s(a[(c0+1)], c[(-2*N-c0+c4-2)]));
          }
          for( c9 = max(-3 * N - c0 + c4 - 1, -16 * c5 - 14); c9 <= min(N, -2 * N - c0 + c4 - 2); c9 += 1) {
            for( c13 = max(1, c0 + 16 * c5 + 1); c13 <= min(c0 + 1, c0 + 16 * c5 + 16); c13 += 1)
              m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] ,H[(c0+1)-c13][c9][(-2*N-c0+c4-c9-1)] - 2*W[c13]);
            for( c13 = max(1, 16 * c5 + c9); c13 <= min(c9, 16 * c5 + c9 + 15); c13 += 1)
              m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-c13][(-2*N-c0+c4-c9-1)] - 2*W[c13]);
            for( c13 = max(1, -2 * N - c0 + c4 + 16 * c5 - c9 - 1); c13 <= min(-2 * N - c0 + c4 - c9 - 1, -2 * N - c0 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9][(-2*N-c0+c4-c9-1)-c13] - 2*W[c13]);
            for( c13 = max(1, 16 * c5 + c9); c13 <= min(min(c0 + 1, c9), 16 * c5 + c9 + 15); c13 += 1)
              m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-c13][c9-c13][(-2*N-c0+c4-c9-1)] - W[c13] + s(a[(c0+1)], b[c9]));
            for( c13 = max(1, -2 * N - c0 + c4 + 16 * c5 - c9 - 1); c13 <= min(min(-2 * N - c0 + c4 - c9 - 1, c9), -2 * N - c0 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-c13][(-2*N-c0+c4-c9-1)-c13] - W[c13] + s(b[c9], c[(-2*N-c0+c4-c9-1)]));
            for( c13 = max(1, -2 * N - c0 + c4 + 16 * c5 - c9 - 1); c13 <= min(min(c0 + 1, -2 * N - c0 + c4 - c9 - 1), -2 * N - c0 + c4 + 16 * c5 - c9 + 14); c13 += 1)
              m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-c13][c9][(-2*N-c0+c4-c9-1)-c13] - W[c13] + s(a[(c0+1)], c[(-2*N-c0+c4-c9-1)]));
          }
        }
        if (c1 == c0 && c2 == 0 && c4 >= 16 && c4 >= 3 * N + c0 + 2 && 2 * N + c0 + 16 >= c4) {
          for( c9 = -3 * N - c0 + c4 - 1; c9 <= N; c9 += 1) {
            m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] ,H[(c0+1)-(c0+1)][c9][(-2*N-c0+c4-c9-1)] - 2*W[(c0+1)]);
            m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-c9][(-2*N-c0+c4-c9-1)] - 2*W[c9]);
            m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - 2*W[(-2*N-c0+c4-c9-1)]);
            if (c0 + 1 >= c9)
              m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-c9][c9-c9][(-2*N-c0+c4-c9-1)] - W[c9] + s(a[(c0+1)], b[c9]));
            if (2 * N + c0 + 2 * c9 + 1 >= c4)
              m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-(-2*N-c0+c4-c9-1)][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - W[(-2*N-c0+c4-c9-1)] + s(b[c9], c[(-2*N-c0+c4-c9-1)]));
            if (2 * N + 2 * c0 + c9 + 2 >= c4)
              m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-(-2*N-c0+c4-c9-1)][c9][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - W[(-2*N-c0+c4-c9-1)] + s(a[(c0+1)], c[(-2*N-c0+c4-c9-1)]));
          }
        } else if (c1 == c0 && c2 == 0 && 3 * N + c0 + 1 >= c4 && c0 + 49 >= c4) {
          for( c5 = max(0, (-c0 + c4 - 1) / 16 - 2); c5 <= min(N / 16, c4 / 16 - 1); c5 += 1) {
            if (c5 == 0)
              for( c9 = 1; c9 < -2 * N - c0 + c4 - 1; c9 += 1) {
                m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] ,H[(c0+1)-(c0+1)][c9][(-2*N-c0+c4-c9-1)] - 2*W[(c0+1)]);
                m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-c9][(-2*N-c0+c4-c9-1)] - 2*W[c9]);
                m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - 2*W[(-2*N-c0+c4-c9-1)]);
                if (c0 + 1 >= c9)
                  m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-c9][c9-c9][(-2*N-c0+c4-c9-1)] - W[c9] + s(a[(c0+1)], b[c9]));
                if (2 * N + c0 + 2 * c9 + 1 >= c4)
                  m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-(-2*N-c0+c4-c9-1)][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - W[(-2*N-c0+c4-c9-1)] + s(b[c9], c[(-2*N-c0+c4-c9-1)]));
                if (2 * N + 2 * c0 + c9 + 2 >= c4)
                  m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-(-2*N-c0+c4-c9-1)][c9][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - W[(-2*N-c0+c4-c9-1)] + s(a[(c0+1)], c[(-2*N-c0+c4-c9-1)]));
              }
            for( c9 = max(max(-c0 + (c0 + c4) / 2 - 8, -N - c0 + (N + c0 + c4) / 2), -c0 - 8 * c5 + (c0 + c4 + 1) / 2 - 8); c9 <= min(min(16, N), -c0 - 8 * c5 + (c0 + c4 + 1) / 2 - 1); c9 += 1)
              m5[(c0+1)][c9][(-c0+c4-2*c9-1)] = INT_MIN;
          }
        }
        for( c9 = max(max(16 * c2 + 1, -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 17), -3 * N - c1 + c4 - 1); c9 <= min(min(N, 16 * c2 + 16), -2 * N - 16 * c0 + 15 * c1 + 16 * c2 + c4 - 2); c9 += 1) {
          if (c1 == c0 && c2 == 0 && c4 <= 15) {
            m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m1[(c0+1)][c9][(-2*N-c0+c4-c9-1)] ,H[(c0+1)-(c0+1)][c9][(-2*N-c0+c4-c9-1)] - 2*W[(c0+1)]);
            m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m2[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-c9][(-2*N-c0+c4-c9-1)] - 2*W[c9]);
            m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m3[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - 2*W[(-2*N-c0+c4-c9-1)]);
            if (c0 + 1 >= c9)
              m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m4[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-c9][c9-c9][(-2*N-c0+c4-c9-1)] - W[c9] + s(a[(c0+1)], b[c9]));
            if (2 * N + c0 + 2 * c9 + 1 >= c4)
              m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m5[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)][c9-(-2*N-c0+c4-c9-1)][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - W[(-2*N-c0+c4-c9-1)] + s(b[c9], c[(-2*N-c0+c4-c9-1)]));
            if (2 * N + 2 * c0 + c9 + 2 >= c4)
              m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)] = MAX(m6[(c0+1)][c9][(-2*N-c0+c4-c9-1)], H[(c0+1)-(-2*N-c0+c4-c9-1)][c9][(-2*N-c0+c4-c9-1)-(-2*N-c0+c4-c9-1)] - W[(-2*N-c0+c4-c9-1)] + s(a[(c0+1)], c[(-2*N-c0+c4-c9-1)]));
          }
          H[(c1+1)][c9][(-2*N-c1+c4-c9-1)] = MAX(0, MAX( H[(c1+1)-1][c9-1][(-2*N-c1+c4-c9-1)-1] + s(a[(c1+1)], b[c9]) + s(a[(c1+1)], c[(-2*N-c1+c4-c9-1)]) + s(b[c9], c[(-2*N-c1+c4-c9-1)]), MAX(m1[(c1+1)][c9][(-2*N-c1+c4-c9-1)], MAX(m2[(c1+1)][c9][(-2*N-c1+c4-c9-1)], MAX(m3[(c1+1)][c9][(-2*N-c1+c4-c9-1)], MAX(m4[(c1+1)][c9][(-2*N-c1+c4-c9-1)], MAX(m5[(c1+1)][c9][(-2*N-c1+c4-c9-1)], m6[(c1+1)][c9][(-2*N-c1+c4-c9-1)])))))));
        }
        if (c1 == c0 && c2 == 0 && c4 <= 15)
          for( c9 = -N - c0 + (N + c0 + c4) / 2; c9 <= N; c9 += 1)
            m5[(c0+1)][c9][(-c0+c4-2*c9-1)] = INT_MIN;
      }
    }
// ---------------------------------------




}