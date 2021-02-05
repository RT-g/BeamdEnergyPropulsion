# BEPのプログラム共有用のリポジトリ
研究成果を統合しやすくしたり、解析に使いたいコードがあちこちに行かないように管理する。

# コードの変数の命名規則
必ず変数は以下のルールに従うこと。
- m_e = 9.10938356e-31 Electron mass
- q = 1.60217662e-19;                        % Elementary charge
- e = 8.85418782e-12;                        % Permittivity
- u = 4.0*3.1415e-7;                         % Permeability
- c = (1.0/(e*u))^0.5;                       % Speed of ligh
- kB = 1.380649e-23;                         % Boltzmann's constant
- Avogadro = 6.02e+23;                       % Avogadro's constant
- h = 6.62607015e-34;                        % Planck's constant
- m_O2 = 32/Avogadro*1e-3;                   % N2 mass
- m_N2 = 28/Avogadro*1e-3;                   % O2 mass
- ratio_N2 = 0.8;                            % Mass ratio of N2
- ratio_O2 = 0.2;                            % Mass ratio of O2
- m_air = m_N2 * ratio_N2 + m_O2 * ratio_O2;     % Air mass
- R_air = kB/m_air;                          % Gas constant
- e_vib_N2 = 0.288;        % 窒素分子の振動エネルギー[eV]<br>(2016_Electron-Vibrational energy exchange), 3380 K(実在気体流れの力学)
- e_vib_O2; % 酸素分子の振動エネルギー[eV]<br>
(2016_Electron-Vibration relaxation in oxygen plasma), 2230 K(実在気体流れの力学)
- e_D; 窒素分子の乖離エネルギー[ev]

## Milliemter-wave and Plasma
- f = 170e+9;                                % Frequency
- omega = 2 * pi * f;                            % Angular frequency
- lambda = c/f;                              % Wavelength