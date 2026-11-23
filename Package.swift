// swift-tools-version: 6.4
import PackageDescription

let package = Package(
    name: "swift-eigen",
    platforms: [.macOS(.v27), .iOS(.v27), .tvOS(.v27), .watchOS(.v27), .visionOS(.v27)],
    products: [
        .library(name: "Eigen", targets: ["Eigen"]),
        .library(name: "Eigen Test Support", targets: ["Eigen Test Support"]),
    ],
    dependencies: [
        .package(url: "https://github.com/swift-atoms/swift-cardinal.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-sparse.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-threefry.git", branch: "main"),
    ],
    targets: [
        .target(
            name: "Eigen",
            dependencies: [
                .product(name: "Cardinal", package: "swift-cardinal"),
                .product(name: "Sparse", package: "swift-sparse"),
                .product(name: "Threefry", package: "swift-threefry"),
            ]
        ),
        .target(
            name: "Eigen Test Support",
            dependencies: [
                .target(name: "Eigen"),
                .product(name: "Sparse Test Support", package: "swift-sparse"),
            ]
        ),
        .testTarget(
            name: "Eigen Tests",
            dependencies: [
                .target(name: "Eigen"),
                .target(name: "Eigen Test Support"),
                .product(name: "Sparse Test Support", package: "swift-sparse"),
            ]
        ),
    ],
    swiftLanguageModes: [.v6]
)

for target in package.targets where ![.system, .binary, .plugin, .macro].contains(target.type) {
    target.swiftSettings =
        (target.swiftSettings ?? []) + [
            .strictMemorySafety(),
            .enableUpcomingFeature("ExistentialAny"),
            .enableUpcomingFeature("InternalImportsByDefault"),
            .enableUpcomingFeature("MemberImportVisibility"),
            .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
            .enableUpcomingFeature("InferIsolatedConformances"),
            .enableExperimentalFeature("Lifetimes"),
            .treatAllWarnings(as: .error),
        ]
}
