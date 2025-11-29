clear all;
clc;
close all;

%% --- GENEL AYARLAR ---
% Bu n değerini değiştirdiğinde hem Kinematik hem Dinamik kısım güncellenir.
n = 2; 
fprintf('=== %d SERBESTLİK DERECELİ (DOF) ROBOT İÇİN HESAPLAMA BAŞLIYOR ===\n', n);

%% --- BÖLÜM 1: İLERİ KİNEMATİK (DH Yöntemi ile Konum Kontrolü) ---
fprintf('\n--- BÖLÜM 1: DH Matrisleri ile Uç Nokta Konumları ---\n');

% Sembolik Değişkenler (Kinematik İçin)
a_dh = sym('L', [1 n]);       % Link uzunlukları (L1, L2...)
theta_dh = sym('theta', [1 n]); % Eklem açıları (theta1, theta2...)

% DH Parametreleri (Planar robot için alpha ve d genelde 0'dır)
alpha = zeros(1, n); 
d = zeros(1, n);     

% Kümülatif Dönüşüm Matrisi (Başlangıçta Birim Matris)
T_cumulative = sym(eye(4)); 

% Döngü ile Her Linkin Hesaplanması
for i = 1:n
    % 1. Mevcut eklemin Dönüşüm Matrisini (A_i) hesapla
    ct = cos(theta_dh(i));
    st = sin(theta_dh(i));
    ca = cos(alpha(i));
    sa = sin(alpha(i));
    
    % DH Matrisi (Standart)
    A_i = [ ct    -st*ca   st*sa   a_dh(i)*ct;
            st     ct*ca  -ct*sa   a_dh(i)*st;
            0      sa      ca      d(i);
            0      0       0       1 ];
        
    % 2. Kümülatif Matrisi Güncelle (T_0i = T_0(i-1) * A_i)
    T_cumulative = T_cumulative * A_i;
    
    % Sadeleştirme
    T_cumulative = simplify(T_cumulative, 'Steps', 10);
    
    % 3. O anki linkin uç noktasını (x, y) çek
    x_pos = T_cumulative(1, 4);
    y_pos = T_cumulative(2, 4);
    
    % 4. Sonuçları Yazdır
    fprintf('Link %d Uç Noktası (x,y):\n', i);
    fprintf('x%d = %s\n', i, char(x_pos));
    fprintf('y%d = %s\n', i, char(y_pos));
end
disp('--- Kinematik Matris Hesabı Tamamlandı ---');
disp('------------------------------------------');

%% --- BÖLÜM 2: DİNAMİK (Lagrange - Enerji Hesabı) ---
fprintf('\n--- BÖLÜM 2: Lagrange Metodu ile Dinamik Denklemler ---\n');

% Sembolik Değişkenler (Dinamik İçin)
syms g real % Yerçekimi
q = sym('q', [n 1], 'real');   % Eklem Pozisyonları (theta yerine q kullanıldı)
dq = sym('dq', [n 1], 'real'); % Eklem Hızları (theta_dot)
m = sym('m', [n 1], 'real');   % Kütleler
L = sym('L', [n 1], 'real');   % Uzunluklar (L1, L2...)

% Enerji Toplamlarını Sıfırla
T_total = 0; % Toplam Kinetik Enerji
P_total = 0; % Toplam Potansiyel Enerji

% Kinematik Takibi İçin Değişkenler
x_prev = 0;      % Bir önceki eklemin x konumu (Taban 0)
y_prev = 0;      % Bir önceki eklemin y konumu (Taban 0)
theta_cum = 0;   % Kümülatif (Toplam) Açı
omega_cum = 0;   % Kümülatif Açısal Hız

% Döngü (Her Uzuv İçin Enerji Hesabı)
for i = 1:n
    % --- KİNEMATİK ---
    % Açıyı güncelle (Absolute Angle) - q(i) burada relative açıdır
    theta_cum = theta_cum + q(i);
    
    % Açısal hızı güncelle (önceki hızlar + şimdiki hız)
    omega_cum = omega_cum + dq(i);
    
    % Kütle Merkezi (CoM) Konumu (Düzgün çubuk: L/2)
    % Not: Önceki eklem konumu + (L/2 * açı)
    x_cm = x_prev + (L(i)) * cos(theta_cum);
    y_cm = y_prev + (L(i)) * sin(theta_cum);
    
    % Bir sonraki eklem için "Eklem Noktası"nı güncelle (Tam boy L kadar git)
    x_prev = x_prev + L(i) * cos(theta_cum);
    y_prev = y_prev + L(i) * sin(theta_cum);
    
    % --- HIZ HESABI (Jacobian Mantığı) ---
    % Kütle merkezinin x ve y hızlarını (türevi) buluruz: v = J * dq
    % Konumun q'ya göre türevi * dq bize hızı verir (Zincir Kuralı)
    vx_cm = jacobian(x_cm, q) * dq; 
    vy_cm = jacobian(y_cm, q) * dq;
    
    % Hızın Karesi (v^2 = vx^2 + vy^2)
    v_squared = simplify(vx_cm^2 + vy_cm^2);
    
    % --- DİNAMİK PARAMETRELER ---
    % 1. Atalet Momenti (Kütle Merkezi Etrafında, Düzgün Çubuk: 1/12)
    % Eğer noktasal kütle istersen burayı 0 yapabilirsin.
    %I_cm =  1*m(i) * L(i)^2;
     I_cm =  0;
    % 2. Kinetik Enerji (Öteleme + Dönme)
    % T = 1/2 * m * v_cm^2  +  1/2 * I_cm * omega_cm^2
    T_link = 0.5 * m(i) * v_squared + 0.5 * I_cm * omega_cum^2;
    T_total = T_total + T_link;
    
    % 3. Potansiyel Enerji (P = m * g * h_cm)
    % Not: Dikey düzlemde çalıştığımızı varsayıyoruz (y ekseni yükseklik)
    P_link = m(i) * g * y_cm;
    P_total = P_total + P_link;
end

%% 3. Sonuçları Göster
fprintf('\n--- SONUÇLAR ---\n');
fprintf('Toplam Kinetik Enerji (T):\n');
% 'collect' komutu denklemi m(1), m(2)... terimlerine göre toparlar
T_grouped = collect(simplify(T_total), m);
P_grouped = collect(simplify(P_total), m);
pretty(T_grouped);

fprintf('\nToplam Potansiyel Enerji (P) [m parantezinde]:\n');
pretty(P_grouped);

fprintf('\nLagrangian (L = T - P):\n');
Lagrangian = simplify(T_total - P_total);
% pretty(Lagrangian); % Çok uzun olabilir, istersen yorumdan çıkar
disp('Lagrangian hesaplandı. (Değişken adı: Lagrangian)');

%% --- BÖLÜM 3: EULER-LAGRANGE İLE TORK HESABI ---
fprintf('\n--- BÖLÜM 3: Tork Denklemlerinin Türetilmesi ---\n');

% İvme değişkenlerini tanımlayalım (ddq: theta_double_dot)
ddq = sym('ddq', [n 1], 'real'); 
tau = sym('tau', [n 1], 'real'); % Tork vektörü

for i = 1:n
    fprintf('Hesaplanıyor: Tork %d...\n', i);
    
    % ADIM 1: Lagrangian'ın HIZA göre kısmi türevi (dL / d_dq)
    dL_d_dqi = diff(Lagrangian, dq(i));
    
    % ADIM 2: Zamana göre türev alma (d/dt)
    % Not: MATLAB sembolik değişkende "t" olmadığı için Zincir Kuralı uygularız:
    % dF/dt = (dF/dq * dq) + (dF/ddq * ddq)
    term_time_derivative = jacobian(dL_d_dqi, q) * dq + jacobian(dL_d_dqi, dq) * ddq;
    
    % ADIM 3: Lagrangian'ın KONUMA göre kısmi türevi (dL / d_q)
    dL_d_qi = diff(Lagrangian, q(i));
    
    % ADIM 4: Tork Denklemi (Euler-Lagrange)
    tau(i) = term_time_derivative - dL_d_qi;
    
    % Sonucu sadeleştir
    tau(i) = simplify(tau(i));
end

%% --- SONUÇLARI YAZDIR ---
fprintf('\n==============================================\n');
fprintf('   GENEL DİNAMİK TORK DENKLEMLERİ (tau)   \n');
fprintf('==============================================\n');
for i = 1:n
    fprintf('\n>> Eklem %d Tork Denklemi (tau_%d) = \n', i, i);
    pretty(tau(i));
end

%% --- BÖLÜM 4: TORK MATRİSİNİN (VEKTÖRÜNÜN) YAZDIRILMASI ---
fprintf('\n==============================================\n');
fprintf('        GENEL TORK VEKTÖRÜ [tau]              \n');
fprintf('==============================================\n');
% Tork vektörünü (tau) ekranda göster
Torque_Matrix = tau; 
fprintf('--- Matematiksel Görünüm ---\n');
pretty(Torque_Matrix);


%% --- BÖLÜM 5: M, C (veya N), G MATRİSLERİNİN ÇIKARILMASI ---
% Denklem Formu: M(q)*ddq + C(q,dq)*dq + G(q) = tau
% Bazen N(q,dq) = C*dq + G olarak gösterilir.

fprintf('\n=======================================================\n');
fprintf('   DİNAMİK MATRİSLERİN AYRIŞTIRILMASI (M, N, G)        \n');
fprintf('   Standart Form: M(q)*ddq + N(q,dq) = tau             \n');
fprintf('   Burada N = C_vektoru + G_vektoru şeklindedir.       \n');
fprintf('=======================================================\n');

% 1. G (Gravity) Vektörü:
% Hız (dq) ve İvme (ddq) sıfır olursa, denklemde sadece Yerçekimi kalır.
G_matrix = subs(Torque_Matrix, [dq; ddq], [zeros(n,1); zeros(n,1)]);
G_matrix = simplify(G_matrix);

% 2. M (Inertia/Kütle) Matrisi:
% Tork denkleminin ivmeye (ddq) göre türevi (Jacobian) bize M matrisini verir.
M_matrix = jacobian(Torque_Matrix, ddq);
M_matrix = simplify(M_matrix);

% 3. C_vektör (Coriolis ve Merkezkaç Terimleri):
% Toplam Torktan, (M*ddq) ve G'yi çıkarırsak geriye hız terimleri kalır.
C_vector = simplify(Torque_Matrix - (M_matrix * ddq) - G_matrix);

% 4. N Matrisi (Lineer Olmayan Terimler):
% Genellikle N = C_vektör + G_vektör olarak istenir.
N_matrix = simplify(C_vector + G_matrix);

% --- SONUÇLARI GÖSTER ---

fprintf('\n>>> M (Kütle/Atalet) Matrisi [%dx%d] <<<\n', n, n);
pretty(M_matrix);

fprintf('\n>>> G (Yerçekimi) Vektörü [%dx1] <<<\n', n);
pretty(G_matrix);

fprintf('\n>>> C (Coriolis & Merkezkaç) Vektörü [%dx1] <<<\n', n);
% Bu vektör C(q,dq)*dq çarpımının sonucudur.
pretty(C_vector);

fprintf('\n>>> N (Non-Linear) Matrisi/Vektörü [%dx1] <<<\n', n);
% N = C_vektör + G
pretty(N_matrix);


%% --- BÖLÜM 5 & 6: DİNAMİK MATRİSLERİN HESAPLANMASI VE SIMULINK ÇIKTISI ---
fprintf('\n==========================================================\n');
fprintf('   1. MATRİSLER HESAPLANIYOR (M, N, G) ...\n');
% ---------------------------------------------------------
% Önceki adımlardan 'tau', 'ddq', 'dq' ve 'n' değişkenlerinin 
% workspace'te tanımlı olduğundan emin olun.
% ---------------------------------------------------------

% 1. M(q) - Atalet Matrisi Hesabı
% Tork denkleminin ivmelere (ddq) göre türevi M matrisini verir.
M_mat = jacobian(tau, ddq);
M_mat = simplify(M_mat);

% 2. G(q) - Yerçekimi Vektörü Hesabı
% Hızları (dq) ve İvmeleri (ddq) sıfır kabul edersek geriye G kalır.
G_vec = subs(tau, [dq; ddq], [zeros(n,1); zeros(n,1)]);
G_vec = simplify(G_vec);

% 3. N(q,dq) - Coriolis ve Merkezkaç Vektörü Hesabı
% Toplam Torktan, (M * ivme) ve (G) çıkarılırsa geriye N kalır.
N_vec = simplify(tau - (M_mat * ddq) - G_vec);

fprintf('   2. SIMULINK KODU OLUŞTURULUYOR ...\n');
fprintf('==========================================================\n');
fprintf('%% --- AŞAĞIDAKİ KODU SIMULINK MATLAB FUNCTION İÇİNE YAPIŞTIRIN ---\n\n');

% --- GİRİŞ DEĞİŞKENLERİNİ EŞLEŞTİRME ---
% Sembolik hesaplama q1, q2 ürettiği için, Simulink vektör girişini (q) bunlara atamalıyız.
fprintf('%% -- Giriş Vektörlerini Skalerlere Ata --\n');
for i = 1:n
    fprintf('q%d = q(%d);\n', i, i);
    fprintf('dq%d = dq(%d);\n', i, i);
end
fprintf('\n');

% --- M MATRİSİ ÇIKTISI ---
fprintf('%% --- M (Atalet) Matrisi ---\n');
fprintf('M = zeros(%d, %d);\n', n, n);
for r = 1:n
    for c = 1:n
        eqn_str = string(M_mat(r,c));
        fprintf('M(%d,%d) = %s;\n', r, c, eqn_str);
    end
end

% --- N VEKTÖRÜ ÇIKTISI ---
fprintf('\n%% --- N (Coriolis/Merkezkaç) Vektörü ---\n');
fprintf('N = zeros(%d, 1);\n', n);
for r = 1:n
    eqn_str = string(N_vec(r));
    fprintf('N(%d) = %s;\n', r, eqn_str);
end

% --- G VEKTÖRÜ ÇIKTISI ---
fprintf('\n%% --- G (Yerçekimi) Vektörü ---\n');
fprintf('G = zeros(%d, 1);\n', n);
for r = 1:n
    eqn_str = string(G_vec(r));
    fprintf('G(%d) = %s;\n', r, eqn_str);
end

fprintf('\n%% --- ÇIKIŞ TORKU (Computed Torque Control Yasası Örneği) ---\n');
fprintf('%% tau = M * (ddq_ref + Kp*e + Kd*de) + N + G;\n');
fprintf('==========================================================\n');





