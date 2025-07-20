/*
    This file is part of OpenTSS.
    Copyright (C) 2022 LatticeX Foundation.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
use crate::FE;
use classgroup::gmp::mpz::Mpz;
use classgroup::gmp_classgroup::*;
use classgroup::ClassGroup;
use curv::arithmetic::Converter;
use curv::arithmetic::*;
use curv::elliptic::curves::Scalar;
use curv::BigInt;
use lazy_static::lazy_static;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CLGroup {
    // fundamental discriminant with $$\Delta_k \equiv 1 \pmod{4}$$.
    // rust analyzer发现该字段只被赋值, 没有被读取.
    pub delta_k: Mpz,

    pub generator: GmpClassGroup,
    pub stilde: Mpz,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PK(pub GmpClassGroup);

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Ciphertext {
    pub c1: GmpClassGroup,
    pub c2: GmpClassGroup,
}

impl From<PK> for GmpClassGroup {
    fn from(pk: PK) -> Self {
        pk.0
    }
}

impl From<GmpClassGroup> for PK {
    fn from(cl: GmpClassGroup) -> Self {
        Self(cl)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SK(pub Mpz);

impl From<SK> for Mpz {
    fn from(sk: SK) -> Self {
        sk.0
    }
}

impl From<Mpz> for SK {
    fn from(mpz: Mpz) -> Self {
        Self(mpz)
    }
}

lazy_static! {
    pub static ref DISCRIMINANT_1827: Mpz = Mpz::from_str("-75257495770792601579408435348799912112609846029965206820064851604692987230254538914853608976971793980958712372789231634579578971529235823075608739231635687425758158575368321348137900869894119507551586698602273331769113654968615517566745786072923103207661147676790644792111452136974276225728730910712947503901232735129687891775293591232029998265064837518833536297518857716272011348573253397254136847763813364524813537416619588617528698171849359403663703760169261184343946919401092992684996593982744033815507830560787451354075275532210193117085590501285653650352846925182015277946751628767130269342252523310043345421861896214174850131607385236887381965429994384214519104490505249675175386383257705274311668138257554180057201072703457873180274207162029503126883077609392094864657038777406276133886450239").unwrap();
}

///According to the paper https://eprint.iacr.org/2020/196.pdf, Sutherland's group-order algorithm may have an impact on security level (1827-bit discriminant cannot ensure 128-bit security). Sutherland's algorithm can compute the order of group and the runtime depends on the order itself rather than the supposed
///size of the order, and hence there is a set of weak instances. The class group generated randomly is vulnerable with small probability(at least 2^{-14} in group with 1827-bit discriminant) and cannot be checked because the order is known, but as the discriminant increases, the probability of successful attack will decrease.
/// The table below shows the relationship between the security level (lambda, rho) and the group size.
/// |lambda\rho|40     |55     |64     |80     |100   |128  |
/// |55                  |660   |825   |880  |1045|1265|1430|
/// |80                  |960   |1200|1280|1520|1840|2080|
/// |100                |1200|1500|1600|1900|2300|2600|
/// |128                |1536|1920|2048|2432|2944|3392|
/// security level (lambda,rho) means the attacker succeeds with probability 1/2^{rho} requires 2^{lambda} bit operations.
/// We provides two groups options: 1827-bit discriminant and 3072-bit discriminant.
impl CLGroup {
    pub fn new_1827() -> Self {
        // 是[CL15, Proposition 1]中的 $$\Delta_k$$.
        // 已验证过是两个质数 $$p, q$$ 的乘积.
        // 参数选取的原理详见[CL15, Appendix B.3].
        let delta_k = Mpz::from_str("-5612960460354297586496608465355436736175385121665162536528003724349027131555226649274328061478036486426974235182817460231858406454328229705097433539599357659030732986212902896965288623752937699627896244889952350312271535460213196686033784826094098560791044370859682930856242386198578254852455887200105136848768296981731378965699234956909793269449142655809687632817484368532297652832818925682445449730939672558315001010323704348812542103398759340104715127787089082447127193712577594846384285770469931817870736146192486488946997648500323172668328291265422577316785106221217309556660122713505680384876843920057653776862871100907889289236674725514431").unwrap();

        let a = Mpz::from_str("3379933361837959750444281267886081834476751587152191195702130129876229099797314884670653751744957540137083102210369145718831424083421213040698452363387299065826090566614550509104171596193940708452801446727936908797340323098201338663853170233065328696856790082422069275092967399794413723895514088363951458374936750806184395472544267780653575123461655052057240595359404437943529185106860238910043016082").unwrap();
        let b = Mpz::from_str("58358596530709071629230628954813789065094567413901151732504604054459961302465715041370372364950254062052414177175583619344532154277172761099891464143583046235404103174114873829883081661462607082144282568946995469931366172071928031362252538721358169137643386731728896321136677327778862260030176007687015790858390775199286445826383171957023481318023285705914617463624817890014105071550499557399120835").unwrap();

        // 是[CL15, Proposition 1]中的 $$\Delta_p=\Delta_k p^2$$.
        let discriminant = Mpz::from_str("-75257495770792601579408435348799912112609846029965206820064851604692987230254538914853608976971793980958712372789231634579578971529235823075608739231635687425758158575368321348137900869894119507551586698602273331769113654968615517566745786072923103207661147676790644792111452136974276225728730910712947503901232735129687891775293591232029998265064837518833536297518857716272011348573253397254136847763813364524813537416619588617528698171849359403663703760169261184343946919401092992684996593982744033815507830560787451354075275532210193117085590501285653650352846925182015277946751628767130269342252523310043345421861896214174850131607385236887381965429994384214519104490505249675175386383257705274311668138257554180057201072703457873180274207162029503126883077609392094864657038777406276133886450239").unwrap();

        // 生成元是一个二次型的等价类 $$[(a, b, c)]$$, 满足 $$\Delta = b^2-4ac$$.
        // [CL15, Fig. 2] $$g:=\left[ \varphi_p^{-1}(\mathfrak{r}^2) \right]^p f^k$$,
        // 其中 $$k$$ 是 1 到 p-1 之间的任意整数.
        // [CL15, Appendix B.1] 
        // $$\varphi_p: I(\mathcal{O}_{\Delta_p}, p)->I(\mathcal{O}_{\Delta_K}, p), \mathfrak{a} \mapsto \mathfrak{a}\mathcal{O}_{\Delta_K}$$.
        // [CL15, Appendix B.4] $$g$$ 的阶约为类数, 类数约为 $$\sqrt{\Delta_K}$$.
        // 之所以说大约, 是因为计算准确的阶数首先要计算类数, 还要计算质因数分解.
        let generator = ClassGroup::from_ab_discriminant(a, b, discriminant);

        let stilde = Mpz::from_str("70874029964003222178994413383062782755071292199599732976843764646488791400299245173357367622414689715904677764175683692699088623752022377648358556868028456505343659927114861398173913787770528036913753917714784290366762147149325499950491790497996441006302782823370615596812470224184985789821376325103006605987671787325355230432").unwrap();
        Self {
            delta_k,
            generator,
            stilde,
        }
    }

    pub fn new_3072() -> Self {
        // 是[CL15, Proposition 1]中的 $$\Delta_k$$.
        let delta_k = Mpz::from_str("-4059187479482350050615628258855828167626431824732199036597668525464616895922000411261718516567731632732286800934600249406393974357768444047141581621951155803795734117021495676831593033172450357785597776576612281305223919414836213766354372816990863555296830253123574199460146205334642841425167146191511265843560519935132843345652241452096808325636749679870044168299284188041110855817763388520168386219623910310164928704787483081634387756726626535065281682599731277374016734081858737636466840542887162979503417512544889504232167650829937659939952944676065304893114687576168003023224828141758525768773373824222139881461335520424806873120226629820060875152488085708505799289587546695067879685280385374856021956449469249646800629229020371797593504643496190406594392765693007499422572180546825466666141075563827212225011483631613617098804995744522667871405671831585120704467080787250858292339350012462220525281878018038188111302643").unwrap();

        let a = Mpz::from_str("7558696258536269836685598691067254392455432126512045048703153470216587589262584772773889350803142696082644492391440222657387505436226106559387636399090692870409664685458565139268115361587965012801626486214178296724409957196458581564654675171891414734873480981011022434578765292108902615901895824430417440790854617643067937627692652614825761175005451229674634875483251243329253319146864770567796872059102546431266905942844241327140097892962334716147988086934626736364915485808267753549075272883526450647112751694864522351665428653263635553").unwrap();
        let b = Mpz::from_str("7110582194089178292873623928339940221177068678147419045016231011875502636852334272115320612748035025965493444044134858046528138916161220583195510694897710613486866021680645839297391773510359545101415795781046778268789176573129050721030057929435162739737320039997979549217218832394942012916585101410393854693994253094796506466642394248189174331033987416154803874258142847897908542918940158825420984527661422824424966184192353102243120918368701972881794424452449979886520228687483571195987688449901751348213142154870047359909042087991612521").unwrap();

        // 是[CL15, Proposition 1]中的 $$\Delta_p=\Delta_k p^2$$.
        let discriminant = Mpz::from_str("-54424806076527156168985032165443480786420749773635324654757897883565645703417884567871042148029631370681593998050904556738585685752650398063780963909714157075453959236964108559638106174791817697784880543469483206812928940067399311317488294362527938882175697366019957696115809536326918268075743444713094760405458243480351677889162968503597995433089866471083629557434694747572996742811168571602020576064530026860732056884277663707523145711611191348850530092216164179043365835314045243342514552963960414585323924973601164865544638680028963602014477110061372790986953196941681814148922608605077598475440517010589628638951818530996940003106846018822916313104592745851208891211437536842831920118837146579414872086574070852830341252313194989178993890549487609254953085788716425539786601425248044456208912081099634485955683519682483772537685513745962700726289912780998237173705138945626073085145984550221831801591311255789856639888250191333084021087921178381947976571438316855053762214058567794937388453551692847435930047709140590212761974735590043211787370187396675877349868884982570867").unwrap();

        // 生成元是一个二次型的等价类 $$[(a, b, c)]$$, 满足 $$\Delta = b^2-4ac$$.
        // TODO: gq的阶是多少?
        let gene = ClassGroup::from_ab_discriminant(a, b, discriminant);

        // $$\tilde{s}$$ 是CL群的阶的上界.
        let stilde = Mpz::from_str("2731990876498942190907198793351360821924936450827254526077205732808204356440122049083260923320622633917729210455296914797563479318897501367228802354913238349385665287912286300665455086668936692955454575005791947875391212727463655396061046670508369983948246816429384317848036518361689362084276319232647078502064602526624505278574375502609123069687358026449142841870608209630753106164304656955565967571069523057451219185931205895267929394930386842181744618272983847612").unwrap();
        Self {
            delta_k,
            generator: gene,
            stilde,
        }
    }

    // 2025.07.16. 此时的generator是 $$f=(p^2, p)$$ 吗?
    pub fn update_class_group_by_p(group: &CLGroup) -> CLGroup {
        let q: Mpz = q();
        let mut gq_new = group.generator.clone();
        gq_new.pow(q);
        CLGroup {
            delta_k: group.delta_k.clone(),
            generator: gq_new,
            stilde: group.stilde.clone(),
        }
    }

    // 源码 `keygen.rs` 用的是 `GROUP_1827`
    pub fn keygen(&self) -> (SK, PK) {
        let sk = SK(bigint_to_mpz(&BigInt::sample_below(
            &(&(mpz_to_bigint(&self.stilde)) * BigInt::from(2u32).pow(40)),
        )));
        let mut generator = self.generator.clone();
        generator.pow(sk.clone().0);
        let pk = PK(generator);
        (sk, pk)
    }

    // 在源码 `sign.rs` 中, `group` 是 `GROUP_UPDATE_1827`
    pub fn encrypt(group: &CLGroup, public_key: &PK, m: &FE) -> (Ciphertext, SK) {
        let m = into_mpz(m);
        let (r, r_big) = group.keygen();
        let delta = group.generator.discriminant().clone();
        let exp_f = expo_f(&q(), &delta, &m);
        let mut h_exp_r = public_key.0.clone();
        h_exp_r.pow(r.0.clone());

        // [CL15, Fig. 1] $$h=g^x, c_1=g^r, c_2=f^mh^r$$.
        let ct = Ciphertext {
            c1: r_big.0,
            c2: h_exp_r * exp_f,
        };
        (ct, r)
    }

    pub fn decrypt(group: &CLGroup, secret_key: &SK, c: &Ciphertext) -> FE {
        // $$(c_1^x)^{-1} == g^{-xr} == h^{-r}$$.
        let mut c1_x_inv = c.c1.clone();
        c1_x_inv.pow(secret_key.0.clone());
        c1_x_inv.inverse();

        // 用 `c1_x_inv` 消掉 $$h^r$$.
        let tmp = c.c2.clone() * &c1_x_inv;

        // 调用离散对数函数, 解出明文.
        let plaintext = discrete_log_f(&q(), &group.generator.discriminant(), &tmp);
        debug_assert!(plaintext < q());
        let plaintext_big = BigInt::from_str_radix(&plaintext.to_str_radix(16), 16).unwrap();
        Scalar::from(&plaintext_big)
    }

    pub fn encrypt_without_r(group: &CLGroup, m: &FE) -> (Ciphertext, SK) {
        let r = SK::from(Mpz::from(0));
        let r_big = group.pk_for_sk(r.clone());
        let m_mpz = Mpz::from_str(&m.to_bigint().to_str_radix(10)).unwrap();
        let exp_f = expo_f(&q(), &group.generator.discriminant(), &m_mpz);

        (
            Ciphertext {
                c1: r_big.0,
                c2: exp_f,
            },
            r,
        )
    }

    pub fn pk_for_sk(&self, sk: SK) -> PK {
        let mut group_element = self.generator.clone();
        group_element.pow(sk.0);
        PK(group_element)
    }

    pub fn eval_scal(c: &Ciphertext, val: Mpz) -> Ciphertext {
        let mut c1 = c.c1.clone();
        c1.pow(val.clone());
        let mut c2 = c.c2.clone();
        c2.pow(val);
        let c_new = Ciphertext { c1, c2 };
        c_new
    }

    pub fn eval_sum(c1: &Ciphertext, c2: &Ciphertext) -> Ciphertext {
        let c_new = Ciphertext {
            c1: c1.c1.clone() * c2.c1.clone(),
            c2: c1.c2.clone() * c2.c2.clone(),
        };
        c_new
    }
}

// secp256k1曲线群的阶
pub fn q() -> Mpz {
    let q = Mpz::from_str(&FE::group_order().to_str_radix(10)).unwrap();
    q
}

// 根据 [Binary Quadratic Forms An Algorithmic Approach, Definition 5.3.3],
// 该函数是在生成 (指定判别式的) 理想类群的零元.
// 零元作为一个理想类, 类内的理想都是主理想 (单一生成理想), 类外的理想都不是主理想.
// [TODO: 2025.07.16. 找到出处证明上述逻辑关系].
pub fn principal_ideal_class(delta: &Mpz) -> GmpClassGroup {
    assert_eq!(delta.mod_floor(&Mpz::from(4)), Mpz::one());
    assert!(delta < &Mpz::zero()); // in general delta can be positive but we don't deal with that case
    let a = Mpz::one();
    let b = Mpz::one();
    ClassGroup::from_ab_discriminant(a, b, (*delta).clone())
}

// 令 $$f=(p^2, p)$$, 计算 $$f^m$$ 的化简形式.
// 本实现并没有先表示出 $$f$$, 再计算 $$f^k$$, 最后化简它;
// 而是利用 [CL15, Proposition 1] 中的公式
// $$\mathtt{Red}(f^k)=(p^2, L(k)p)$$.
pub fn expo_f(p: &Mpz, delta: &Mpz, k: &Mpz) -> GmpClassGroup {
    if k == &Mpz::zero() {
        let group = principal_ideal_class(delta);
        return group;
    }
    let mut k_inv = k.invert(p).unwrap();
    if k_inv.mod_floor(&Mpz::from(2)) == Mpz::zero() {
        k_inv = k_inv - p;
    };
    let k_inv_p = k_inv * p;

    let qf = ClassGroup::from_ab_discriminant(p * p, k_inv_p, (*delta).clone());
    qf
}

pub fn discrete_log_f(p: &Mpz, delta: &Mpz, fm: &GmpClassGroup) -> Mpz {
    let principal_qf = principal_ideal_class(delta);
    if fm == &principal_qf {
        return Mpz::zero();
    } else {
        // `lk` 就是 [CL15, Proposition 1] 中的 $$L(m)$$,
        // 同时 `lk` 又是二次型的 b 参数.
        let lk = fm.b.div_floor(p);
        let lk_inv = lk.invert(p).unwrap();
        return lk_inv;
    }
}

pub fn mpz_to_bigint(value: &Mpz) -> BigInt {
    BigInt::from_str_radix(&value.to_str_radix(16), 16).unwrap()
}

pub fn bigint_to_mpz(value: &BigInt) -> Mpz {
    Mpz::from_str_radix(&value.to_str_radix(16), 16).unwrap()
}

pub fn into_mpz(f: &FE) -> Mpz {
    Mpz::from_str(&f.to_bigint().to_str_radix(10)).unwrap()
}

lazy_static! {
    // [CL15, Fig. 2]
    // $$g:=\left[ \varphi_p^{-1}(\mathfrak{r}^2) \right]^p f^k$$.
    // $$G:=\left< g \right>$$
    pub static ref GROUP_1827: CLGroup = CLGroup::new_1827();
}

lazy_static! {
    // [CL15, Fig. 2]
    pub static ref GROUP_UPDATE_1827: CLGroup = CLGroup::update_class_group_by_p(&GROUP_1827);
}

lazy_static! {
    pub static ref GROUP_3072: CLGroup = CLGroup::new_3072();
}

lazy_static! {
    pub static ref GROUP_UPDATE_3072: CLGroup = CLGroup::update_class_group_by_p(&GROUP_3072);
}

// #[test]
// pub fn test_expo_f() {
//     use curv::elliptic::curves::traits::ECScalar;
//     use curv::BigInt;
//     use curv::arithmetic::Converter;
//     use class_group::BinaryQF;
//     let p_bigint = FE::group_order();
//     let disc_bigint = BigInt::from_str_radix("-427591883024055094237166135622616655692519934789141165516706756107713228024295574699944370011492039773032284231197067920516476468952667833450787437565142878414690959601489749046040711938128287296782189616537717688270214244402769076279798454569992482605974226864130564559661953774942648236131809818008262290621847888715639292230937796387198730350695926159265589312774352775617969025314432860667259678434940630540543035480609747636331455190603326444795763211972135959921012248669035956143099101074949371262559470141902477520013849988014387063025822245215125266145907356818350851633878657454104378492479085550260787918948512444782208939491", 10).unwrap();
//     let k_bigint = BigInt::from_str_radix("12345", 10).unwrap();
//     let result_1 = BinaryQF::expo_f(&p_bigint, &disc_bigint, &k_bigint);
//     println!("result_1 = {:?}", result_1);

//     let p_mpz = Mpz::from_str("115792089237316195423570985008687907852837564279074904382605163141518161494337").unwrap();
//     let disc_mpz = Mpz::from_str("-427591883024055094237166135622616655692519934789141165516706756107713228024295574699944370011492039773032284231197067920516476468952667833450787437565142878414690959601489749046040711938128287296782189616537717688270214244402769076279798454569992482605974226864130564559661953774942648236131809818008262290621847888715639292230937796387198730350695926159265589312774352775617969025314432860667259678434940630540543035480609747636331455190603326444795763211972135959921012248669035956143099101074949371262559470141902477520013849988014387063025822245215125266145907356818350851633878657454104378492479085550260787918948512444782208939491").unwrap();
//     let k_mpz = Mpz::from_str("12345").unwrap();
//     let result_2 = expo_f(&p_mpz, &disc_mpz, &k_mpz);
//     println!("result_2 = {:?}", result_2);
// }

// #[test]
// pub fn test_compose() {
//     use crate::utilities::class::GROUP_128;
//     use curv::BigInt;
//     use curv::arithmetic::Converter;
//     let a = Mpz::from_str("3379933361837959750444281267886081834476751587152191195702130129876229099797314884670653751744957540137083102210369145718831424083421213040698452363387299065826090566614550509104171596193940708452801446727936908797340323098201338663853170233065328696856790082422069275092967399794413723895514088363951458374936750806184395472544267780653575123461655052057240595359404437943529185106860238910043016082").unwrap();
//     let b = Mpz::from_str("58358596530709071629230628954813789065094567413901151732504604054459961302465715041370372364950254062052414177175583619344532154277172761099891464143583046235404103174114873829883081661462607082144282568946995469931366172071928031362252538721358169137643386731728896321136677327778862260030176007687015790858390775199286445826383171957023481318023285705914617463624817890014105071550499557399120835").unwrap();
//     let discriminant = Mpz::from_str("-75257495770792601579408435348799912112609846029965206820064851604692987230254538914853608976971793980958712372789231634579578971529235823075608739231635687425758158575368321348137900869894119507551586698602273331769113654968615517566745786072923103207661147676790644792111452136974276225728730910712947503901232735129687891775293591232029998265064837518833536297518857716272011348573253397254136847763813364524813537416619588617528698171849359403663703760169261184343946919401092992684996593982744033815507830560787451354075275532210193117085590501285653650352846925182015277946751628767130269342252523310043345421861896214174850131607385236887381965429994384214519104490505249675175386383257705274311668138257554180057201072703457873180274207162029503126883077609392094864657038777406276133886450239").unwrap();
//     let gq: GmpClassGroup = ClassGroup::from_ab_discriminant(a.clone(), b.clone(), discriminant.clone());
//     let mut d = gq.clone();
//     d.pow(Mpz::from_str("123").unwrap());
//     let mul = gq*d;
//     println!("mul = {:?}", mul);

//     let gq_1 = GROUP_128.gq.clone();
//     let d_1 = GROUP_128.gq.clone();
//     let d_2 = d_1.exp(&BigInt::from_str_radix("123", 10).unwrap());
//     let comp = gq_1.compose(&d_2).reduce();
//     println!("comp = {:?}", comp);
// }

#[test]
pub fn test_encrypt_decrypt() {
    let start_1827 = time::now();
    let m = FE::random();
    let (sk, pk) = GROUP_1827.keygen();
    let c = CLGroup::encrypt(&GROUP_1827, &pk, &m);
    let m_new = CLGroup::decrypt(&GROUP_1827, &sk, &c.0);
    assert_eq!(m, m_new);
    let end_1827 = time::now();
    println!("time with 1827bit = {:?}", end_1827 - start_1827);
    let start_3072 = time::now();
    let m = FE::random();
    let (sk, pk) = GROUP_3072.keygen();
    let c = CLGroup::encrypt(&GROUP_3072, &pk, &m);
    let m_new = CLGroup::decrypt(&GROUP_3072, &sk, &c.0);
    assert_eq!(m, m_new);
    let end_3072 = time::now();
    println!("time with 3072bit = {:?}", end_3072 - start_3072);
}

#[test]
pub fn pow_a() {
    use crate::GE;
    use anyhow::format_err;
    let mut a = GROUP_1827.generator.clone();
    println!("a = {}", a.to_bytes().len());
    let b = Mpz::from_str_radix("123", 10).unwrap();
    a.pow(b);
    println!("a_new = {}", a.to_bytes().len());
    let m = FE::random();
    let n = bincode::serialize(&m)
        .map_err(|why| format_err!("bincode serialize error: {}", why))
        .unwrap();
    println!("m = {}", m.to_bigint().to_bytes().len());
    println!("n = {}", n.len());
    let a = BigInt::from_bytes(&vec![21; 33]);
    let e: FE = Scalar::from(&a);
    println!("e = {:?}", e);

    let z = e.to_bigint().to_bytes();
    println!("e = {}", e.to_bigint().to_bytes().len());
    let f: FE = Scalar::from(&BigInt::from_bytes(&z));
    println!("f = {:?}", f);

    let s = GE::zero();
    println!("s ={:?}", s);
    let g = s.to_bytes(true);
    let s_1 = GE::from_bytes(&g[1..g.len()]);
    println!("s_1 ={:?}", s_1);

    let c = bincode::serialize(&s)
        .map_err(|why| format_err!("bincode serialize error: {}", why))
        .unwrap();
    println!("c ={}", c.len());
}

#[test]
fn test_big_to_mpz() {
    let a = BigInt::from_str_radix("123", 16).unwrap();
    let start = time::now();
    let _b = bigint_to_mpz(&a);
    let end = time::now();
    println!("duration = {:?}", end - start);
}
