# Mini Plasma Physics（迷你等离子体物理）

**从零开始、零依赖的 C 语言实现**，涵盖大学级等离子体物理与受控核聚变。每个模块对应 MIT、Princeton、Stanford、Cambridge 及其他顶尖大学的一门或多门课程，覆盖从基础等离子体理论（德拜屏蔽、动理学理论、磁流体力学）到聚变能源、空间物理、激光等离子体相互作用、工业等离子体加工、尘埃（复杂）等离子体以及波与不稳定性理论的全谱系。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-plasma-fundamental](mini-plasma-fundamental/) | 德拜屏蔽、等离子体频率/参数、动理学理论（Vlasov–Fokker–Planck）、单/双流体 MHD、粒子模拟（Boris 推进器、PIC）、等离子体波色散 | MIT 22.611, Princeton PHY 521, Berkeley PHYS 242 |
| [mini-mhd](mini-mhd/) | 理想/电阻 MHD 方程、MHD 平衡（Grad–Shafranov、Z-箍缩、θ-箍缩）、MHD 不稳定性（扭曲、腊肠、交换、气球模）、MHD 波模（Alfvén、快/慢磁声波）、MHD 数值方法（Riemann 求解器、约束传输、∇·B 清洁） | MIT 22.611/22.615, Princeton AST 554, Cambridge Part III MHD |
| [mini-waves-instabilities](mini-waves-instabilities/) | 色散关系求解器（Brent、复牛顿、Nyquist）、动理学色散（Landau 阻尼、回旋阻尼、Bernstein 模）、非线性波（三波耦合、Zakharov 方程、孤子）、不稳定性目录（双流、Buneman、离子声、漂移波） | MIT 22.611, Princeton PHY 521, Berkeley PHYS 242 |
| [mini-dusty-plasma](mini-dusty-plasma/) | 尘埃颗粒充电（OML 理论）、尘埃受力（离子拖曳、中性拖曳、热泳、Lorentz 力）、尘埃晶体形成（Coulomb 耦合参数、相变）、尘埃等离子体波模（DAW、DIAW、DLW）、尘埃输运系数 | MIT 22.611/22.612, Princeton PHY 535 |
| [mini-fusion-plasma](mini-fusion-plasma/) | 磁约束（托卡马克几何、磁面、Shafranov 位移）、MHD 平衡/稳定性（Grad–Shafranov、撕裂模、破裂）、加热与电流驱动（NBI、ICRF、ECRH、LH）、聚变输运（新经典、湍流）、Lawson 判据与功率平衡、ITER/DEMO/SPARC 运行点 | MIT 22.611/22.615, Princeton PHY 563, Cambridge Part III Fusion |
| [mini-space-plasma](mini-space-plasma/) | 行星磁层（Chapman–Ferraro、Dungey 重联）、太阳风（Parker 模型、日球层结构）、空间 MHD、空间等离子体波（whistler、chorus、ULF）、空间天气参数 | MIT 22.611/22.616, Cambridge Part III Plasma Physics |
| [mini-laser-plasma](mini-laser-plasma/) | 激光吸收机制（逆轫致辐射、共振吸收、Brunel 加热）、参量不稳定性（SRS、SBS、TPD、丝化）、场致/碰撞电离（BSI、ADK、雪崩）、相对论自聚焦、尾场加速（LWFA）、激光场中单粒子运动（Boris 积分器） | MIT 22.611, Princeton PHY 525, Stanford PHYSICS 370 |
| [mini-industrial-plasma](mini-industrial-plasma/) | RF 放电模型（CCP/ICP）、电子能量分布函数（Maxwell、Druyvesteyn）、等离子体鞘层（Bohm 判据、Child–Langmuir）、等离子体-表面相互作用（刻蚀、PECVD、溅射）、碰撞截面与速率系数 | MIT 22.611, Stanford EE 414, Berkeley EECS 245 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明与参考推导
- **实用演示程序** — Grad–Shafranov 平衡求解器、CMA 图、粒子推进器、鞘层模型、波色散根、不稳定性增长率计算器等

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-plasma-fundamental
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-plasma-physics/
├── mini-plasma-fundamental/    # 基础等离子体物理（参数、动理学、MHD、粒子、波）
├── mini-mhd/                   # 磁流体力学（方程、平衡、不稳定性、数值方法、波）
├── mini-waves-instabilities/   # 等离子体波与不稳定性（色散求解器、动理学、非线性、目录）
├── mini-dusty-plasma/          # 尘埃（复杂）等离子体物理（充电、受力、晶体、波、输运）
├── mini-fusion-plasma/         # 聚变等离子体物理（约束、平衡、加热、MHD、输运）
├── mini-space-plasma/          # 空间等离子体物理（磁层、太阳风、MHD、波）
├── mini-laser-plasma/          # 激光等离子体相互作用（吸收、不稳定性、电离、尾场）
└── mini-industrial-plasma/     # 工业等离子体加工（放电、EEDF、鞘层、表面）

```

## 许可证

MIT
